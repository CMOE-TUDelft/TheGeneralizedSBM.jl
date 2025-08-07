"""
StokesParams

Parameters for the linear stokes problem.
"""
@with_kw struct StokesParams
  domain::Tuple{Float64,Float64,Float64,Float64} = (0.0,1.0,0.0,1.0)
  n_cells::Tuple{Int64,Int64} = (10,10)
  in_out_wall_tags::Tuple = ([7],[8],[1,2,3,4,5,6])
  ϕ::Function = level_set(CircleParams())
  order::Int64 = 1 # FE order
  weight_quad_degree::Int64 = 20 # Quadrature degree for computing weights
  weight_min_value::Float64 = 0.0 # Minimum value for weights
  weight_approach::Symbol = :standard # Approach for computing weights (:standard or :binary (0.0 or 1.0))
  output_folder::String = datadir("sims","Stokes") # Output folder
  λ::Float64 = 0.5 # Fraction for the "weight_approach" :fraction case
  verbose::Bool = true # Verbosity
  f::Function = x -> VectorValue(0.0,0.0) # Source term
  g::Function = n -> VectorValue(0.0,0.0) # Neumann boundary condition
  u₀::Function = x -> VectorValue(0.0,0.0) # Exact solution velocity
  p₀::Function = x -> 0.0 # Exact solution pressure
  ν::Float64 = 1.0 # viscosity
  β₁::Float64 = 20.0 # Penalty parameter
  βᵤ::Float64 = 1.0 # Penalty parameter
  βₚ::Float64 = 0.1 # Penalty parameter
  Cinv::Vector{Float64} = [36.0,36.0] # Inverse estimate constant per order
end

"""
    main_stokes()

Main function for the stokes problem. It computes the solution to the stokes problem
on an arbitrary shaped domain with Dirichlet boundary conditions. The solution is computed using
the TheGeneralizedSBM method, assuming the boundary is defined through a level set function.
"""
function main_stokes(params::StokesParams)

  # General parameters
  @unpack output_folder, verbose = params

  # Discrete model
  @unpack domain, n_cells, in_out_wall_tags = params
  model = CartesianDiscreteModel(domain,n_cells)
  Ω = Interior(model)
  Λ = Skeleton(Ω)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"inlet",in_out_wall_tags[1])
  add_tag_from_tags!(labels,"outlet",in_out_wall_tags[2])
  add_tag_from_tags!(labels,"wall",in_out_wall_tags[3])
  Γ_out = Boundary(Ω,tags="outlet")

  # Compute weights
  α = compute_weights(Ω,params)
  params_ref = reconstruct(params,weight_approach=:binary)
  α_ref = compute_weights(Ω,params_ref)

  # Define FE spaces
  @unpack order = params
  @unpack f, g, u₀, p₀, ν = params
  D = num_dims(Ω)
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)#,space=:P)
  V = TestFESpace(Ω,reffeᵤ,conformity=:H1,dirichlet_tags=["inlet","wall"])
  Q = TestFESpace(Ω,reffeₚ,conformity=:C0)#,constraint=:zeromean)
  U = TrialFESpace(V,[u₀,u₀])
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V,Q])
  X = MultiFieldFESpace([U,P])

  # Define measures
  dΩ = Measure(Ω,2*order)
  dΛ = Measure(Λ,2*order)
  nΛ = get_normal_vector(Λ)
  dΓ = Measure(Γ_out,2*order)
  nΓ = get_normal_vector(Γ_out)

  # Auxiliar variables
  Id = TensorValue(Matrix(1.0I,D,D))
  σ(ε) = 2ν*ε
  h = mean(CellField(lazy_map(vol->(vol)^(1/D),get_cell_measure(Ω)),Ω))
  h3 = mean(CellField(lazy_map(vol->((vol)^(1/D))^3,get_cell_measure(Ω)),Ω))
  αₚ = 1.0 * (mean(α) .> 0.0)  * (mean(α) .< 1.0) * (abs(jump(α)) .< 1.0)
  αₑ = 1.0 * (mean(α).== 0.0)
  dcf = CellField(d(params),Ω)
  u₀d = CellField(x->u₀(x+d(params)(x)),Ω)
  sᵤ(u) = u+dcf⋅∇(u)+1/2*(dcf⋅(dcf⋅∇∇(u)))
  sᵤ¹(u) = u+dcf⋅∇(u)
  verbose &&  writevtk(Λ,"weights_lambda",cellfields=["αp"=>αₚ,"αe"=>αₑ])

  # Define variational problem
  @unpack β₁, βᵤ, βₚ, Cinv = params
  η = 1.0#1/(8*√(Cinv[order]))*(-4+√(Cinv[order])+(65*Cinv[order]+56*√(Cinv[order])+16)^(1/2))
  γ = β₁*(order+1)^2*η*ν # rectangles
  # γ = β₁*(order+1)*(order+2)/2*η*ν #triangres
  κ = 1.0
  βdiv = 1.0
  
  a((u,p),(v,q)) = ∫( α*( σ(ε(u))⊙ε(v) ) )dΩ - ∫( α*( p*(∇⋅v) ) )dΩ + ∫( κ*α*( q*(∇⋅u) ) )dΩ +
    ∫( βdiv*( (∇⋅v)*(∇⋅u) ) )dΩ -
    ∫( meanᵧ(α,σ(ε(u)) - p*Id)⊙jump(α*nΛ⊗v) )dΛ -
    ∫( (meanᵧ(α,σ(ε(v)) + κ*q*Id)⋅jump(α*nΛ))⋅meanᵧ(α,sᵤ(u)) )dΛ +
    ∫( ((γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v))))⋅(meanᵧ(α,sᵤ(u))) )dΛ + 
    ∫( βᵤ*(ν)*(αₚ+αₑ)*( h*(jump(nΛ⋅ε(u))⊙jump(nΛ⋅ε(v))) + h3*((0.5*(nΛ.⁺⋅(jump(nΛ⋅∇∇(u))+jump(nΛ⋅∇∇(u))')))⋅(0.5*(nΛ.⁺⋅(jump(nΛ⋅∇∇(v))+jump(nΛ⋅∇∇(v))')))) ) )dΛ +
    ∫( κ*βₚ/(ν+βdiv)*(αₚ+αₑ)*( h3*(jump(nΛ⋅∇(p))⊙jump(nΛ⋅∇(q))) ) )dΛ 
  l((v,q)) =
    ∫(α*(f⋅v))dΩ + ∫( v⋅g )dΓ -
    ∫( (meanᵧ(α,σ(ε(v)) + κ*q*Id)⋅jump(α*nΛ))⋅u₀d )dΛ +
    ∫( ((γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v)⋅u₀d))) )dΛ
  op = AffineFEOperator(a,l,X,Y)

  # solution
  ls = LUSolver()
  uₕ,pₕ = solve(ls,op)

  # Postprocess
  eᵤl₂ = uₕ - u₀
  eₚl₂ = pₕ-p₀
  l2normᵤ = √(∑(∫(α*(eᵤl₂⋅eᵤl₂))dΩ))
  l2normₚ = √(∑(∫(α*(eₚl₂⋅eₚl₂))dΩ))
  l2normᵤ_ref = √(∑(∫(α_ref*(eᵤl₂⋅eᵤl₂))dΩ))
  l2normₚ_ref = √(∑(∫(α_ref*(eₚl₂⋅eₚl₂))dΩ))

  verbose && writevtk(Ω,output_folder*"/solution",cellfields=["weights"=>α,"uₕ"=>uₕ,"u₀"=>u₀,"pₕ"=>pₕ,"p₀"=>p₀,"dcf"=>dcf,"u₀d"=>u₀d,"ϕ"=>params.ϕ])

  return l2normᵤ,l2normₚ,l2normᵤ_ref,l2normₚ_ref

end

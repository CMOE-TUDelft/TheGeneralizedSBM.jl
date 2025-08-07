"""
LinearElasticityParams

Parameters for the linear elasticity problem.

  * `domain`: Tuple with the domain limits `(xmin,xmax,ymin,ymax)`.
  * `n_cells`: Tuple with the number of cells `(nx,ny)`.
  * `ϕ`: Level set function defining the boundary.
  * `order`: Order of the Finite Element space.
  * `weight_quad_degree`: Quadrature degree for computing weights.
  * `weight_min_value`: Minimum value for weights.
  * `weight_approach`: Approach for computing weights (`:standard`,  `:binary` (0.0 or 1.0) or `:fraction` (optimal surrogate)).
  * `λ`: Fraction for the "weight_approach" :fraction case.
  * `output_folder`: Output folder.
  * `verbose`: Verbosity.
  * `f`: Source term.
  * `u₀`: Exact solution.
  * `μ`: Deviatoric Lame parameter.
  * `λₑ`: Volumetric Lame parameter.
  * `β₁`: Penalty parameter.
  * `β₂`: Penalty parameter.
  * `β₃`: Penalty parameter.
"""
@with_kw struct LinearElasticityParams
  domain::Tuple{Float64,Float64,Float64,Float64} = (0.0,1.0,0.0,1.0)
  n_cells::Tuple{Int64,Int64} = (10,10)
  ϕ::Function = level_set(CircleParams())
  order::Int64 = 1 # FE order
  weight_quad_degree::Int64 = 20 # Quadrature degree for computing weights
  weight_min_value::Float64 = 0.0 # Minimum value for weights
  weight_approach::Symbol = :standard # Approach for computing weights (:standard or :binary (0.0 or 1.0))
  λ::Float64 = 0.5 # Fraction for the "weight_approach" :fraction case
  output_folder::String = datadir("sims","LinearElasticity") # Output folder
  verbose::Bool = true # Verbosity
  f::Function = x -> VectorValue(0.0,0.0) # Source term
  u₀::Function = x -> VectorValue(0.0,0.0) # Exact solution
  μ::Float64 = 1.0 # Deviatoric Lame parameter
  λₑ::Float64 = 1.0 # Volumetric Lame parameter
  β₁::Float64 = 20.0 # Penalty parameter
  β₂::Float64 = 0.5 # Penalty parameter
  β₃::Float64 = 0.5 # Penalty parameter
end

"""
    main_elasticity()

Main function for the linear elasticity problem. It computes the solution to the linear elasticity problem
on an arbitrary shaped domain with Dirichlet boundary conditions. The solution is computed using
the TheGeneralizedSBM method, assuming the boundary is defined through a level set function.
"""
function main_elasticity(params::LinearElasticityParams)

  # General parameters
  @unpack output_folder, verbose = params

  # Discrete module
  @unpack domain, n_cells = params
  # model = simplexify(CartesianDiscreteModel(domain,n_cells))
  model = CartesianDiscreteModel(domain,n_cells)
  Ω = Interior(model)
  Λ = Skeleton(Ω)

  # Compute weights
  α = compute_weights(Ω,params)
  params_ref = reconstruct(params,weight_approach=:binary)
  α_ref = compute_weights(Ω,params_ref)

  # Define FE spaces
  @unpack order = params
  D = num_dims(Ω)
  reffe = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  V = TestFESpace(Ω,reffe,conformity=:H1)
  U = TrialFESpace(V)

  # Define measures
  dΩ = Measure(Ω,2*order)
  dΛ = Measure(Λ,2*order)
  nΛ = get_normal_vector(Λ)

  # Auxiliar variables
  @unpack f, u₀, μ, λₑ = params
  Id = TensorValue(Matrix(1.0I,D,D))
  σ(ε) = 2μ*ε + λₑ*tr(ε)*Id
  D = num_dims(Ω)
  h = mean(CellField(lazy_map(vol->(vol)^(1/D),get_cell_measure(Ω)),Ω))
  h3 = mean(CellField(lazy_map(vol->((vol)^(1/D))^3,get_cell_measure(Ω)),Ω))
  # αₚ = 1/(λₑ+2μ) * (mean(α) .< 1.0) * ( abs(jump(α)) .!= 2*mean(α) )
  # αₑ = 1/(λₑ+2μ) * (mean(α).== jump(α))
  αₚ = 1.0 * (mean(α) .> 0.0)  * (mean(α) .< 1.0) * (abs(jump(α)) .< 1.0)# (α.⁺ .< 1.0) * (α.⁻ .< 1.0) 
  αₑ = 1.0 * (mean(α).== 0.0)
  dcf = CellField(d(params),Ω)
  u₀d = CellField(x->u₀(x+d(params)(x)),Ω)
  sᵤ(u) = u+dcf⋅∇(u)+1/2*(dcf⋅(dcf⋅∇∇(u)))
  sᵤ¹(u) = u+dcf⋅∇(u)

  # Define variational problem
  @unpack β₁, β₂, β₃ = params
  Cinv = 36.0
  η = 1.0#1/(8*√(Cinv))*(-4+√(Cinv)+(65*Cinv+56*√(Cinv)+16)^(1/2))
  γ = β₁*(order+1)^2*η*(λₑ+2μ) #rectangles
  # γ = β₁*(order+1)*(order+2)/2*η*(λₑ+2μ) #triangres
  a(u,v) =
    ∫(α*(σ(ε(u))⊙ε(v)))dΩ -
    ∫( meanᵧ(α,σ(ε(u)))⊙jump(α*v⊗nΛ) +
       meanᵧ(α,σ(ε(v)))⊙(jump(α*nΛ)⊗meanᵧ(α,sᵤ(u))) -
       (γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v))⋅meanᵧ(α,sᵤ(u))) -
       (β₂*αₚ+β₃*αₑ)*(λₑ+2μ)*h*(jump(nΛ⋅∇(u))⊙jump(nΛ⋅∇(v)) ) -
       (β₂*αₚ+β₃*αₑ)*(λₑ+2μ)*h3*( (nΛ.⁺⋅(jump(nΛ⋅∇∇(u)))) ⋅ (nΛ.⁺⋅(jump(nΛ⋅∇∇(v)))) ))dΛ
      #  (β₂*αₚ+β₃*αₑ)*h*(jump(∇(u))⊙jump(∇(v)) ) -
      #  (β₂*αₚ+β₃*αₑ)*h3*( ((jump(∇∇(u)))) ⊙ ((jump(∇∇(v)))) ))dΛ
  l(v) =
    ∫(α*(f⋅v))dΩ -
    ∫( meanᵧ(α,σ(ε(v)))⊙(jump(α*nΛ)⊗meanᵧ(α,u₀d)) -
      (γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v))⋅meanᵧ(α,u₀d)) )dΛ
  op = AffineFEOperator(a,l,V,U)

  # solution
  uₕ = solve(op)
  verbose && writevtk(Ω,output_folder*"/solution",cellfields=["weights"=>α,"uₕ"=>uₕ,"u₀"=>u₀,"dcf"=>dcf,"ϕ"=>params.ϕ])

  # Postprocess
  el₂ = uₕ - u₀
  l2norm = √(∑(∫(α*(el₂⋅el₂))dΩ))

  return l2norm

end

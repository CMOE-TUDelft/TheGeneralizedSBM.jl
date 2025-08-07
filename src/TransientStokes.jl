"""
TransientStokesParams

Parameters for the linear stokes problem.
"""
@with_kw struct TransientStokesParams
  domain::Tuple{Float64,Float64,Float64,Float64} = (0.0,1.0,0.0,1.0)
  n_cells::Tuple{Int64,Int64} = (10,10)
  ϕ::Function = t->level_set(CircleParams())
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
  β₁::Float64 = 1.0 # Penalty parameter
  βᵤ::Float64 = 1.0 # Penalty parameter
  βₚ::Float64 = 0.1 # Penalty parameter
  βdiv::Float64 = 0.0 # Grad-div stabilization parameter
  global_gp::Bool = false # Global ghost penalty stabilization
  Cinv::Vector{Float64} = [36.0,36.0] # Inverse estimate constant per order
  ode_solver_params::TimeIntegratorParams = TimeIntegratorParams()
  is_exact_solution::Bool = false # Use exact solution as initial condition
  in_out_wall_tags::Tuple = ([7],[8],[1,2,3,4,5,6])
  dirichlet_tags::Vector{String} = ["inlet","wall"]
  dirichlet_masks::Vector{Tuple} = [(true,true),(true,true)]
  dirichlet_functions::Vector{Function} = [u₀,u₀]
end

"""
    params(t)

Return the parameters for the transient stokes problem at time t.
"""
(params::TransientStokesParams)(t) = reconstruct(params,ϕ=params.ϕ(t))

"""
    main_transient_stokes()

Main function for the transient stokes problem. It computes the solution to the stokes problem
on an arbitrary shaped domain with Dirichlet boundary conditions. The solution is computed using
the TheGeneralizedSBM method, assuming the boundary is defined through a time-dependent level set function.
"""
function main_transient_stokes(params::TransientStokesParams)

  # General parameters
  @unpack output_folder, verbose = params

  # Discrete model
  @unpack domain, n_cells, in_out_wall_tags = params
  model = CartesianDiscreteModel(domain,n_cells)
  Ω = Interior(model)
  Λ = Skeleton(Ω)
  labels = get_face_labeling(model)
  length(in_out_wall_tags[1])>0 && add_tag_from_tags!(labels,"inlet",in_out_wall_tags[1])
  length(in_out_wall_tags[2])>0 && add_tag_from_tags!(labels,"outlet",in_out_wall_tags[2])
  length(in_out_wall_tags[3])>0 && add_tag_from_tags!(labels,"wall",in_out_wall_tags[3])
  Γ_out = Boundary(Ω,tags="outlet")

  # Compute weights
  α_buffer = Ref{Any}((α=nothing,t=nothing))
  function α(t)
    if α_buffer[].t == t
      return α_buffer[].α
    else
      println("Computing α at t=$t")
      α_buffer[] = (α=compute_weights(Ω,params(t)),t=t)
      println("Done computing α")
      return α_buffer[].α
    end
  end
  α_ref_buffer = Ref{Any}((α=nothing,t=nothing))
  function α_ref(t)
    if α_ref_buffer[].t == t
      return α_ref_buffer[].α
    else
      params_ref = reconstruct(params(t),weight_approach=:binary)
      α_ref_buffer[] = (α=compute_weights(Ω,params_ref),t=t)
      return α_ref_buffer[].α
    end
  end

  # Define FE spaces
  @unpack order = params
  @unpack f, g, u₀, p₀, ν, dirichlet_tags, dirichlet_functions, dirichlet_masks = params
  D = num_dims(Ω)
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)#,space=:P)
  V = TestFESpace(Ω,reffeᵤ,conformity=:H1,dirichlet_tags=dirichlet_tags,dirichlet_masks=dirichlet_masks)
  Q = TestFESpace(Ω,reffeₚ,conformity=:C0)#,constraint=:zeromean)
  U = TransientTrialFESpace(V,dirichlet_functions)
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V,Q])
  X = TransientMultiFieldFESpace([U,P])

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
  @unpack global_gp = params
  function αₚ(t)
    if !global_gp
      return 1.0 * (mean(α(t)) .> 0.0)  * (mean(α(t)) .< 1.0) * (abs(jump(α(t))) .< 1.0)# (α.⁺ .< 1.0) * (α.⁻ .< 1.0) 
    else
      return 1.0 * (mean(α(t)).!= 0.0)
    end
  end 
  αₑ(t) =  1.0 * (mean(α(t)).== 0.0)
  
  dcf(t) = CellField(d(params(t)),Ω)
  u₀d(t) = CellField(x->u₀(x+d(params(t))(x),t),Ω)
  sᵤ(u,t) = u+dcf(t)⋅∇(u)+1/2*(dcf(t)⋅(dcf(t)⋅∇∇(u)))
  sᵤ¹(u,t) = u+dcf(t)⋅∇(u)
  verbose &&  writevtk(Λ,output_folder*"/weights_lambda",cellfields=["αp"=>αₚ(0.0),"αe"=>αₑ(0.0)])

  # Define variational problem
  @unpack β₁, βᵤ, βₚ, βdiv, Cinv = params
  η = 1.0#1/(8*√(Cinv[order]))*(-4+√(Cinv[order])+(65*Cinv[order]+56*√(Cinv[order])+16)^(1/2))
  γ = β₁*(order+1)^2*η*ν # rectangles
  # γ = β₁*(order+1)*(order+2)/2*η*ν #triangres
  κ = 1.0
  
  m(t,(uₜ,),(v,)) = ∫( α(t)*(uₜ⋅v) )dΩ
  a(t,(u,p),(v,q)) = ∫( α(t)*( σ(ε(u))⊙ε(v) ) )dΩ - ∫( α(t)*( p*(∇⋅v) ) )dΩ + ∫( κ*α(t)*( q*(∇⋅u) ) )dΩ +
    ∫( βdiv*( (∇⋅v)*(∇⋅u) ) )dΩ -
    ∫( meanᵧ(α(t),σ(ε(u)) - p*Id)⊙jump(α(t)*nΛ⊗v) )dΛ -
    ∫( (meanᵧ(α(t),σ(ε(v)) + κ*q*Id)⋅jump(α(t)*nΛ))⋅meanᵧ(α(t),sᵤ(u,t)) )dΛ +
    ∫( ((γ/h*abs(jump(α(t))))*(meanᵧ(α(t),sᵤ¹(v,t))))⋅(meanᵧ(α(t),sᵤ(u,t))) )dΛ + 
    ∫( βᵤ*ν*(αₚ(t)+αₑ(t))*( h*(jump(nΛ⋅ε(u))⊙jump(nΛ⋅ε(v))) + h3*((0.5*(nΛ.⁺⋅(jump(nΛ⋅∇∇(u))+jump(nΛ⋅∇∇(u))')))⋅(0.5*(nΛ.⁺⋅(jump(nΛ⋅∇∇(v))+jump(nΛ⋅∇∇(v))')))) ) )dΛ +
    # ∫( κ*βₚ/(ν+βdiv)*(αₚ(t)+αₑ(t))*( h3*(jump(nΛ⋅∇(p))⊙jump(nΛ⋅∇(q))) ) )dΛ 
    ∫( κ*βₚ/(ν)*(αₚ(t)+αₑ(t))*( h3*(jump(nΛ⋅∇(p))⊙jump(nΛ⋅∇(q))) ) )dΛ 
  l(t,(v,q)) =
    ∫(α(t)*(f(t)⋅v))dΩ + ∫( v⋅g(t) )dΓ -
    ∫( (meanᵧ(α(t),σ(ε(v)) + κ*q*Id)⋅jump(α(t)*nΛ))⋅u₀d(t) )dΛ +
    ∫( ((γ/h*abs(jump(α(t))))*(meanᵧ(α(t),sᵤ¹(v,t)⋅u₀d(t)))) )dΛ

  res(t,(u,p),(v,q)) = a(t,(u,p),(v,q)) - l(t,(v,q))
  jac(t,(u,p),(du,dp),(v,q)) = a(t,(du,dp),(v,q))
  jac_t(t,(u,),(duₜ,),(v,)) = m(t,(duₜ,),(v,))
  a₀((u,p),(v,q)) = a(0.0,(u,p),(v,q))
  l₀((v,q)) = l(0.0,(v,q))

  op₀ = AffineFEOperator(a₀,l₀,X(0.0),Y)
  op = TransientLinearFEOperator((a,m),l,X,Y)

  # solution
  ls = LUSolver()
  @unpack is_exact_solution = params
  if is_exact_solution
    xₕ₀ = interpolate_everywhere([u₀(0.0),p₀(0.0)],X(0.0))
    vₕ₀ = interpolate_everywhere([∂t(u₀)(0.0),0.0],X(0.0))
  else
    xₕ₀ = solve(op₀)
    vₕ₀ = interpolate_everywhere([VectorValue(0.0,0.0),0.0],X(0.0))
  end
  verbose && writevtk(Ω,output_folder*"/solution_0.0",cellfields=["weights"=>α(0),"uₕ"=>xₕ₀[1],"u₀"=>u₀(0),"pₕ"=>xₕ₀[2],"p₀"=>p₀(0),"dcf"=>dcf(0),"u₀d"=>u₀d(0),"ϕ"=>params(0).ϕ])

  @unpack ode_solver_params = params
  @unpack T, Δt = ode_solver_params
  ode_solver = get_ode_solver(ode_solver_params,ls)
  xₕ₀solver = get_initial_condition(ode_solver_params,xₕ₀,vₕ₀)
  xₕₜ = solve(ode_solver,op,0.0,T,xₕ₀solver)

  # Norms
  l2norm(a,t) = √(∑(∫(α(t)*(a⋅a))dΩ)) 
  h1norm(a,t) = √(∑(∫(α(t)*(∇(a)⊙∇(a)))dΩ))
  l2l2ᵤ = 0.0
  l2h1ᵤ = 0.0
  l2l2ₚ = 0.0
  l2h1ₚ = 0.0
  l∞l2ᵤ = 0.0
  l∞h1ᵤ = 0.0
  l∞l2ₚ = 0.0
  l∞h1ₚ = 0.0
  global lastl2ᵤ = 0.0
  global lastl2ₚ = 0.0
  FD = Float64[]
  time = Float64[]

  nd(t) = dcf(t)/(dcf(t)⋅dcf(t)).^(0.5)
  nu(t) = VectorValue(1.0,0.0)#u₀d(t)/(u₀d(t)⋅u₀d(t)).^(0.5)

  # Postprocess  
  cnt = 0
  createpvd(output_folder*"/solution") do pvd
    for (t,(uh,ph)) in xₕₜ
      global lastl2ᵤ = l2norm(uh-u₀(t),t)
      global lastl2ₚ = l2norm(ph-p₀(t),t)
      l2l2ᵤ += √(Δt*(l2norm(uh-u₀(t),t))^2)
      l2h1ᵤ += √(Δt*(h1norm(uh-u₀(t),t))^2)
      l2l2ₚ += √(Δt*(l2norm(ph-p₀(t),t))^2)
      l2h1ₚ += √(Δt*(h1norm(ph-p₀(t),t))^2)
      l∞l2ᵤ = max(l∞l2ᵤ,l2norm(uh-u₀(t),t))
      l∞h1ᵤ = max(l∞h1ᵤ,h1norm(uh-u₀(t),t))
      l∞l2ₚ = max(l∞l2ₚ,l2norm(ph-p₀(t),t))
      l∞h1ₚ = max(l∞h1ₚ,h1norm(ph-p₀(t),t))
      push!(FD,∑(∫(abs(jump(α_ref(t)))*meanᵧ(α_ref(t),sᵤ¹(ph,t)*nΛ)*(nd(t)⋅nu(t))⋅nd(t))dΛ))
      push!(time,t)
      cnt += 1
      if verbose 
        pvd[t] = createvtk(Ω,output_folder*"/solution_$t",cellfields=["weights"=>α(t),"uₕ"=>uh,"u₀"=>u₀(t),"pₕ"=>ph,"p₀"=>p₀(t),"dcf"=>dcf(t),"u₀d"=>u₀d(t),"ϕ"=>params(t).ϕ],order=2)
      end
    end
  end

  return l2l2ᵤ,l2h1ᵤ,l2l2ₚ,l2h1ₚ,l∞l2ᵤ,l∞h1ᵤ,l∞l2ₚ,l∞h1ₚ,lastl2ᵤ,lastl2ₚ,FD,time

end

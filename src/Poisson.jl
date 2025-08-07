"""
    PoissonParams

Parameters for the Poisson problem. The parameters are:

    * `domain`: Tuple with the domain limits `(xmin,xmax,ymin,ymax)`.
    * `n_cells`: Tuple with the number of cells `(nx,ny)`.
    * `ϕ`: Level set function defining the boundary.
    * `order`: Order of the Finite Element space.
    * `weight_quad_degree`: Quadrature degree for computing weights.
    * `weight_min_value`: Minimum value for weights.
    * `weight_approach`: Approach for computing weights (`:standard` or `:binary` (0.0 or 1.0)).
    * `output_folder`: Output folder.
    * `verbose`: Verbosity.
    * `f`: Source term.
    * `u₀`: Exact solution.
"""
@with_kw struct PoissonParams
  domain::Tuple{Vararg{Float64}} = (0.0,1.0,0.0,1.0)
  n_cells::Tuple{Vararg{Int64}} = (10,10)
  ϕ::Function = level_set(CircleParams())
  order::Int64 = 1 # FE order
  weight_quad_degree::Int64 = 20 # Quadrature degree for computing weights
  weight_min_value::Float64 = 0.0 # Minimum value for weights
  weight_approach::Symbol = :standard # Approach for computing weights (:standard or :binary (0.0 or 1.0))
  λ::Float64 = 0.5 # Fraction for the "weight_approach" :fraction case
  output_folder::String = datadir("sims","Poisson") # Output folder
  verbose::Bool = true # Verbosity
  compute_cond::Bool = false # Compute condition number
  f::Function = x -> 0.0 # Source term
  g::Function = n -> 0.0 # Neumann boundary condition
  u₀::Function = x -> 0.0 # Exact solution
  β₁::Float64 = 20.0 # Penalty parameter
  β₂::Float64 = 0.5 # Penalty parameter
  β₃::Float64 = 0.5 # Penalty parameter
end

"""
    main_poisson()

Main function for the Poisson problem. It computes the solution to the Poisson problem
on an arbitrary shaped domain with Dirichlet boundary conditions. The solution is computed using
the TheGeneralizedSBM method, assuming the boundary is defined through a level set function.
"""
function main_poisson(params::PoissonParams)

  # General parameters
  @unpack output_folder, verbose = params

  # Discrete module
  @unpack domain, n_cells = params
  # model = simplexify(CartesianDiscreteModel(domain,n_cells))
  model = CartesianDiscreteModel(domain,n_cells)
  Ω = Interior(model)
  Λ = Skeleton(Ω)
  Γ = Boundary(Ω,tags="boundary")

  # Compute weights
  α = compute_weights(Ω,params)
  params_ref = reconstruct(params,weight_approach=:binary)
  α_ref = compute_weights(Ω,params_ref)

  # Define FE spaces
  @unpack order = params
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffe,conformity=:H1)
  U = TrialFESpace(V)

  # Define measures
  dΩ = Measure(Ω,2*order)
  dΛ = Measure(Λ,2*order)
  dΓ = Measure(Γ,2*order)
  nΛ = get_normal_vector(Λ)
  nΓ = get_normal_vector(Γ)

  # Auxiliar variables
  @unpack f, g, u₀ = params
  D = num_dims(Ω)
  h = mean(CellField(lazy_map(vol->(vol)^(1/D),get_cell_measure(Ω)),Ω))
  h2 = mean(CellField(lazy_map(vol->((vol)^(1/D))^2,get_cell_measure(Ω)),Ω))
  h3 = mean(CellField(lazy_map(vol->((vol)^(1/D))^3,get_cell_measure(Ω)),Ω))
  # αₚ = 1.0 * (mean(α) .< 1.0) * ( abs(jump(α)) .!= 2*mean(α) )
  # αₑ = 1.0 * (mean(α).== jump(α))
  αₚ = 1.0 * (mean(α) .> 0.0)  * (mean(α) .< 1.0) * (abs(jump(α)) .< 1.0)# (α.⁺ .< 1.0) * (α.⁻ .< 1.0) 
  αₑ = 1.0 * (mean(α).== 0.0)
  dcf = CellField(d(params),Ω)
  u₀d = CellField(x->u₀(x+d(params)(x)),Ω)
  sᵤ(u) = u+dcf⋅∇(u)+1/2*(dcf⋅(dcf⋅∇∇(u)))
  sᵤ¹(u) = u+dcf⋅∇(u)

  # Define weak form
  @unpack β₁, β₂, β₃ = params
  Cinv = 36.0
  η = 1.0#1/(8*√(Cinv))*(-4+√(Cinv)+(65*Cinv+56*√(Cinv)+16)^(1/2))
  γ = β₁*(order+1)^2*η #rectangles
  # γ = β₁*(order+1)*(order+2)/2*η #triangres
  a(u,v) =
    ∫(α*(∇(u)⋅∇(v)) )dΩ +
    ∫( -1.0*(meanᵧ(α,∇(u))⋅jump(α*v*nΛ)) -
       (meanᵧ(α,∇(v))⋅jump(α*nΛ))*meanᵧ(α,sᵤ(u)) +
      (γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v))*meanᵧ(α,sᵤ(u))) +
      (β₂*αₚ+β₃*αₑ)*h*(jump(∇(u)⋅nΛ)⋅jump(∇(v)⋅nΛ) ) +
      (β₂*αₚ+β₃*αₑ)*h3*( (nΛ.⁺⋅(jump(nΛ⋅∇∇(u)))) ⋅(nΛ.⁺⋅(jump(nΛ⋅∇∇(v))) )))dΛ +
    ∫( (0.0 .< α .< 1.0)*(-1.0*(∇(u)⋅(α*v*nΓ)) - (∇(v)⋅(α*nΓ))*(sᵤ(u)) + (γ/h*α)*(sᵤ¹(v)*sᵤ(u))) )dΓ
  l(v) =
    ∫(α*(f*v))dΩ + ∫(α*(g(nΓ)*v))dΓ +
    ∫( -1.0*(meanᵧ(α,∇(v))⋅jump(α*nΛ))*meanᵧ(α,u₀d) +
      (γ/h*abs(jump(α)))*(meanᵧ(α,sᵤ¹(v))*meanᵧ(α,u₀d)) )dΛ +
    ∫( (0.0 .< α .< 1.0)*(-1.0*(∇(v)⋅(α*nΓ))*u₀d + (γ/h*α)*(sᵤ¹(v)*u₀d)) )dΓ
  op = AffineFEOperator(a,l,V,U)
  A = get_matrix(op)
  if params.compute_cond
    cond_A = cond(Matrix(A),2)
  else
    cond_A = 0.0
  end

  # solution
  uₕ = solve(op)
  verbose && writevtk(Ω,output_folder*"/solution",cellfields=["weights"=>α,"uₕ"=>uₕ,"u₀"=>u₀,"dcf"=>dcf,"ϕ"=>params.ϕ])

  # Postprocess
  el₂ = uₕ - u₀
  l2norm = √(∑(∫(α*el₂*el₂)dΩ))

  return l2norm, cond_A

end

module ManufacturedTransientStokes_dt
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TheGeneralizedSBM
using Gridap
using TimerOutputs
using Plots
using DataFrames, DataFramesMeta

to = TimerOutput("ManufacturedTransientStokesDtConvergence")

# Define fixed parameters (M. Olshanskii)
domain = (-1.0,2.0,-1.0,1.0)
center(t) = [t,0.0]
radius = 0.5
u₁(x,t) = 2π*x[2]*cos(π*((x[1]-t)^2+x[2]^2))
u₂(x,t) = -2π*(x[1]-t)*cos(π*((x[1]-t)^2+x[2]^2))
p(x,t) = sin(π*((x[1]-t)^2+x[2]^2))-2/π
@show umax = √((2π*2)^2 + (2π)^2)
@show hmin = 2/80
@show Δtmin = hmin/umax
@show T = 20Δtmin

u(x,t) = VectorValue(u₁(x,t),u₂(x,t))
u(t::Real) = x -> u(x,t)
p(t::Real) = x -> p(x,t)
ν = 0.01
I = TensorValue(1.0,0.0,0.0,1.0)
f(t::Real) = x -> ∂t(u)(t)(x) - ν*Δ(u(t))(x) + ∇(p(t))(x)
nout = VectorValue(1.0,0.0)
εn(x,t) = nout⋅ε(u(t))(x)
∇un(x,t) = nout⋅∇(u(t))(x)
g(x,t) = 2*ν*εn(x,t) - p(x,t)*nout
g(t::Real) = x -> g(x,t)
β₁=1.0
βᵤ=1.0
βₚ=1.0e-1
order = 2

# Warm-up parameters
ϕ(t) = level_set(CircleParams(center=center(t),radius=radius,in_out=1))
n_cells = (9,6)
output_folder = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedTransientStokes_dt")
ode_solver_params = TimeIntegratorParams(T=0.1,method=:generalized_alpha,ρ∞=0.0)
params = TransientStokesParams(
  domain=domain,
  ϕ=ϕ,
  f=f,
  g=g,
  u₀=u,
  p₀=p,
  ν=ν,
  n_cells=n_cells,
  output_folder=output_folder,
  order=2,
  ode_solver_params=ode_solver_params,
  is_exact_solution=true,
  verbose=false,
)

# Execute main function (Warm-up)
println(main_transient_stokes(params))

# Test parameters dt convergence
weight_approach = [:standard]
nΔt = [4,2,1,0.5,0.25]
ϕ_name = [:circle]
βdiv = [1.0]
global_gp = [false,true]
all_params = @strdict weight_approach nΔt ϕ_name βdiv global_gp
cases = dict_list(all_params)

# Get level set functions
function get_level_set(ϕ_name::Symbol)
  if ϕ_name == :circle
    return t->level_set(CircleParams(center=center(t),radius=radius,in_out=1))
  elseif ϕ_name == :flower
    return t->level_set(FlowerParams(center=center(t),radius=radius,n=n_petals,in_out=1))
  elseif ϕ_name == :parallelogram
    return t->level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=-1))
  else
    error("Case not recognized")
  end
end

# Execute dt convergence case function
function execute_case(case)
  @unpack weight_approach, nΔt, ϕ_name, βdiv, global_gp = case
  case_name = savename(case,"jld2",allowedtypes=(Real, String, Symbol, Function))
  println("Executing case: ",case_name)

  n = 20
  Δt = nΔt * Δtmin

  # Case parameters
  ϕ = get_level_set(ϕ_name)
  if weight_approach == :fraction
    λ = 0.5
  else 
    λ = 1.0
  end
  ode_solver_params = TimeIntegratorParams(Δt=Δt,T=T,method=:generalized_alpha,ρ∞=0.0)
  verbose = false
  output_folder = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedTransientStokes_dt",replace(case_name,".jld2"=>""))
  if !isdir(output_folder)
    mkpath(output_folder)
  end
  params = TransientStokesParams(
    domain=domain,ϕ=ϕ,
    f=f,g=g,
    u₀=u,p₀=p,
    ν=ν,
    n_cells=(3n/2,n),
    weight_approach=weight_approach,
    verbose=verbose,
    order=order,
    β₁=β₁,βᵤ=βᵤ,βₚ=βₚ,βdiv=βdiv,
    global_gp=global_gp,
    output_folder=output_folder,
    λ=λ,
    ode_solver_params=ode_solver_params,
    is_exact_solution=true,
  )
  println("Parameters: ",params)

  # Execute main function
  results = copy(case)
  @timeit to "main_$(case_name)" results["l2l2ᵤ"], results["l2h1ᵤ"], results["l2l2ₚ"], results["l2h1ₚ"], results["l∞l2ᵤ"], results["l∞h1ᵤ"], results["l∞l2ₚ"], results["l∞h1ₚ"], results["lastl2ᵤ"], results["lastl2ₚ"] = main_transient_stokes(params)
  results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

  return results
end

# Execute dt convergence cases
for case in cases
  path = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedTransientStokes_dt")
  filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
  data, file = produce_or_load(path,case,execute_case;filename=filename)
end

# Get data
all_results = collect_results(datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedTransientStokes_dt"))
plot_geom_cases = [:circle]#,:flower,:parallelogram]
plot_approaches = [:standard]#,:binary,:fraction]
plot_βdiv = [1.0]
plot_global_gp = [false,true]

# Plot results
labels = ["WSBM","SBM","λ-SBM"]
markers = [:circle,:square,:utriangle]
lines = [:solid,:dash]
colors = ["#0072B2", "#E69F00", "#009E73"]
plt_label = "Error norm"
plt_name = ["l2l2ᵤ","l2l2ₚ","l∞l2ᵤ","l∞l2ₚ"]
plt_slope = [0.5,0.5,0.0,1.0]
shift = [3.0e-3,2.0e-2,5.0e-1,2.0e-3]
for (igeom, geom_case) in enumerate(plot_geom_cases)
  plts = []
  plot_ref_line = Bool[]
  for iplot in 1:10
    push!(plts,plot(xlabel="Time step size",ylabel=plt_label,legend=:topright,xaxis=:log,yaxis=:log))
    push!(plot_ref_line,true)
  end
  for (iβdiv,βdiv) in enumerate(plot_βdiv)
    for (iglobal_gp,global_gp) in enumerate(plot_global_gp)
      for (iapproach, approach) in enumerate(plot_approaches)
        results = @linq all_results[all_results.:ϕ_name.==geom_case .&& 
                                    all_results.:weight_approach.==approach .&&
                                    all_results.:βdiv.==βdiv .&&
                                    all_results.:global_gp.==global_gp ,:] |> orderby(:nΔt)
        for iplot in 1:length(plt_name)
          if iplot == 3
            plot!(plts[iplot],results.:nΔt*Δtmin,results[!,plt_name[iplot]],marker=markers[iglobal_gp],ls=lines[iβdiv],color=colors[iglobal_gp],label="global_gp=$(global_gp)",ylims=(5e-2,1e0),xticks = (results.:nΔt*Δtmin, string.(results.:nΔt) .* "Δtmin"))
          else
            plot!(plts[iplot],results.:nΔt*Δtmin,results[!,plt_name[iplot]],marker=markers[iglobal_gp],ls=lines[iβdiv],color=colors[iglobal_gp],label="global_gp=$(global_gp)")
          end
          if plot_ref_line[iplot] && iplot !=3
            plot!(plts[iplot],results.:nΔt*Δtmin,shift[iplot]*(results.:nΔt*Δtmin).^(-plt_slope[iplot]),xticks = (results.:nΔt*Δtmin, string.(results.:nΔt) .* "Δtmin"),ls=:dashdot,color=:black,label=false)
            x_triangle = [2Δtmin, 2.5Δtmin, 2.5Δtmin]
            y_triangle = [shift[iplot]*(x_triangle[1]).^(-plt_slope[iplot]),shift[iplot]*(x_triangle[1]).^(-plt_slope[iplot]),shift[iplot]*(x_triangle[2]).^(-plt_slope[iplot])]
            plot!(plts[iplot],x_triangle, y_triangle, lw = 1, color = :black, label = "")
            annotate!(2.2Δtmin, (1.1*shift[iplot])*(x_triangle[1]).^(-plt_slope[iplot]), text("1", :black, 9, :left))
            annotate!(2.6Δtmin, 0.95*shift[iplot]*(2Δtmin).^(-plt_slope[iplot]), text(string(plt_slope[iplot]), :black, 9, :left))
            plot_ref_line[iplot] = false
          end
        end
      end
    end
  end
  for iplot in 1:length(plt_name)
    savefig(plts[iplot],plotsdir("Journal_paper_GSBM","Stokes","ManufacturedTransientStokes_dt",plt_name[iplot]*"_"*string(ϕ_name[igeom])*".pdf"))
  end
end

end
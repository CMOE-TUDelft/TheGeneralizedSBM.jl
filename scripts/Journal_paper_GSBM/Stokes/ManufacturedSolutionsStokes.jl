module ManufacturedSolutionsStokes
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TheGeneralizedSBM
using Gridap
using TimerOutputs
using Plots
using DataFrames, DataFramesMeta

to = TimerOutput("ManufacturedStokes")

# Define fixed parameters
domain = (-1.0,2.0,-1.0,1.0)
center = [0.0,0.0]
radius = 0.5
n_petals = 5
u₁(x) = 2π*x[2]*cos(π*((x[1])^2+x[2]^2))
u₂(x) = -2π*(x[1])*cos(π*((x[1])^2+x[2]^2))
u(x) = VectorValue(u₁(x),u₂(x))
p(x) = sin(π*((x[1])^2+x[2]^2))-2/π
ν = 0.01
I = TensorValue(1.0,0.0,0.0,1.0)
f(x) = -ν*Δ(u)(x) + ∇(p)(x)
nout = VectorValue(1.0,0.0)
εn(x) = nout⋅ε(u)(x)
g(x) = 2*ν*εn(x) - p(x)*nout

# Warm-up parameters
ϕ = level_set(CircleParams(center=center,radius=radius,in_out=1))
n_cells = (12,6)
output_folder = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedStokes")
params = StokesParams(ϕ=ϕ,f=f,g=g,u₀=u,p₀=p,ν=ν,n_cells=n_cells,output_folder=output_folder,order=2,weight_quad_degree=20)

# Execute main function (Warm-up)
println(main_stokes(params))

# Test parameters
weight_approach = [:standard,:binary,:fraction]
verbose = false
n = [10,20,40,80,160]
ϕ_name = [:circle,:flower,:parallelogram]
order = [2]#,2]
all_params = @strdict weight_approach verbose n ϕ_name order
cases = dict_list(all_params)
β₁=1.0
βᵤ=1.0
βₚ=1.0e-1

# Execute case function
function execute_case(case)
  @unpack weight_approach, verbose, n, ϕ_name, order = case
  case_name = savename(case,"jld2",allowedtypes=(Real, String, Symbol, Function))
  println("Executing case: ",case_name)

  # Case parameters
  if ϕ_name == :circle
    ϕ = level_set(CircleParams(center=center,radius=radius,in_out=1))
  elseif ϕ_name == :flower
    ϕ = level_set(FlowerParams(center=center,radius=radius,n=n_petals,in_out=1))
  elseif ϕ_name == :parallelogram
    ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=-1))
  else
    error("Case not recognized")
  end
  if weight_approach == :fraction
    λ = 0.5
  else 
    λ = 1.0
  end
  params = StokesParams(ϕ=ϕ,f=f,g=g,u₀=u,p₀=p,ν=ν,n_cells=(3n/2,n),weight_approach=weight_approach,
    verbose=verbose,order=order,β₁=β₁,βᵤ=βᵤ,βₚ=βₚ,output_folder=datadir("sims","ManufacturedStokes"),λ=λ)

  # Execute main function
  results = copy(case)
  @timeit to "main_$(case_name)" results["l2normᵤ"], results["l2normₚ"], results["l2normᵤ_ref"], results["l2normₚ_ref"] = main_stokes(params)
  results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

  return results
end

# Execute cases
for case in cases
  path = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedStokes")
  filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
  data, file = produce_or_load(path,case,execute_case;filename=filename)
end

# Get data
all_results = collect_results(datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedStokes"))
plot_geom_cases = [:circle,:flower,:parallelogram]
plot_approaches = [:standard,:binary,:fraction]
plot_orders = [2]#,2]

# Plot results
labels = ["WSBM","SBM","λ-SBM"]
markers = [:circle,:square,:utriangle]
lines = [:solid,:dash]
colors = ["#0072B2", "#E69F00", "#009E73"]
filenames_u = ["l2norm_circle_u.pdf","l2norm_flower_u.pdf","l2norm_parallelogram_u.pdf"]
filenames_p = ["l2norm_circle_p.pdf","l2norm_flower_p.pdf","l2norm_parallelogram_p.pdf"]
filenames_u_ref = ["l2norm_circle_u_ref.pdf","l2norm_flower_u_ref.pdf","l2norm_parallelogram_u_ref.pdf"]
filenames_p_ref = ["l2norm_circle_p_ref.pdf","l2norm_flower_p_ref.pdf","l2norm_parallelogram_p_ref.pdf"]
shift_u = 1.0e2
shift_p = 1.0e0
for (igeom, geom_case) in enumerate(plot_geom_cases)
  plt = plot(xlabel="Number of elements per direction",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
  plt_ref = plot(xlabel="Number of elements per direction",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
  plt2 = plot(xlabel="Number of elements per direction",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
  plt2_ref = plot(xlabel="Number of elements per direction",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
  plot_ref_line = true
  for (iorder, order) in enumerate(plot_orders)
    for (iapproach, approach) in enumerate(plot_approaches)
      results = @linq all_results[all_results.:ϕ_name.==geom_case .&& all_results.:weight_approach.==approach .&& all_results.:order.==order,:] |> orderby(:n)
      plot!(plt,results.:n,results.:l2normᵤ,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
      plot!(plt2,results.:n,results.:l2normₚ,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
      plot!(plt_ref,results.:n,results.:l2normᵤ_ref,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
      plot!(plt2_ref,results.:n,results.:l2normₚ_ref,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
      if plot_ref_line
        plot!(plt,results.:n,shift_u*(results.:n).^(-3),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
        plot!(plt2,results.:n,shift_p*(results.:n).^(-2),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
        plot!(plt_ref,results.:n,shift_u*(results.:n).^(-3),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
        plot!(plt2_ref,results.:n,shift_p*(results.:n).^(-2),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
        x_triangle = [80, 110, 110]
        y_triangle = [shift_u*(x_triangle[1]).^(-3),shift_u*(x_triangle[1]).^(-3),shift_u*(x_triangle[2]).^(-3)]
        y_triangle2 = [shift_p*(x_triangle[1]).^(-2),shift_p*(x_triangle[1]).^(-2),shift_p*(x_triangle[2]).^(-2)]
        plot!(plt,x_triangle, y_triangle, lw = 1, color = :black, label = "")
        annotate!(90, 1.3*shift_u*(x_triangle[1]).^(-3), text("1", :black, 9, :left))
        annotate!(115, shift_u*(90).^(-3), text("3", :black, 9, :left))
        plot!(plt2,x_triangle, y_triangle2, lw = 1, color = :black, label = "")
        annotate!(90, 1.3*shift_p*(x_triangle[1]).^(-2), text("1", :black, 9, :left))
        annotate!(115, shift_p*(90).^(-2), text("2", :black, 9, :left))
        plot!(plt_ref,x_triangle, y_triangle, lw = 1, color = :black, label = "")
        annotate!(90, 1.3*shift_u*(x_triangle[1]).^(-3), text("1", :black, 9, :left))
        annotate!(115, shift_u*(90).^(-3), text("3", :black, 9, :left))
        plot!(plt2_ref,x_triangle, y_triangle2, lw = 1, color = :black, label = "")
        annotate!(90, 1.3*shift_p*(x_triangle[1]).^(-2), text("1", :black, 9, :left))
        annotate!(115, shift_p*(90).^(-2), text("2", :black, 9, :left))
        plot_ref_line = false
      end
    end
  end
  savefig(plt,plotsdir("Journal_paper_GSBM","Stokes","ManufacturedStokes",filenames_u[igeom]))
  savefig(plt2,plotsdir("Journal_paper_GSBM","Stokes","ManufacturedStokes",filenames_p[igeom]))
  savefig(plt_ref,plotsdir("Journal_paper_GSBM","Stokes","ManufacturedStokes",filenames_u_ref[igeom]))
  savefig(plt2_ref,plotsdir("Journal_paper_GSBM","Stokes","ManufacturedStokes",filenames_p_ref[igeom]))
end

function plot_vtk(n)
  n_cells = (3n/2,n)

  # VTK results circle case
  ϕ = level_set(CircleParams(center=center,radius=radius,in_out=1))
  output_folder = datadir("sims","Journal_paper_GSBM","Stokes","ManufacturedStokes","VTK_circle")
  params = StokesParams(domain=domain,ϕ=ϕ,f=f,g=g,u₀=u,p₀=p,ν=ν,order=2,n_cells=n_cells,output_folder=output_folder,β₁=β₁,βᵤ=βᵤ,βₚ=βₚ,weight_approach=:fraction,λ=0.5)
  println(main_stokes(params))

end

plot_vtk(40)

end

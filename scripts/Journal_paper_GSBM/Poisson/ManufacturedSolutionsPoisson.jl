module ManufacturedSolutionsPoisson
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TheGeneralizedSBM
using Gridap
using TimerOutputs
using Plots
using DataFrames, DataFramesMeta

to = TimerOutput("ManufacturedPoisson")

# Define fixed parameters
center = [0.5,0.5]
radius = 0.3
n_petals = 5
u(x) = sin(15*π*x[1])*sin(15*π*x[2])
f(x) = -tr(∇∇(u)(x))

# Warm-up parameters
ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=1))
n_cells = (8,8)
output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson")
params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,order=2)

# Execute main function (Warm-up)
println(main_poisson(params))

# Test parameters
weight_approach = [:standard,:binary,:fraction]
verbose = false
n = [20,40,80,160,320]
order = [1]
ϕ_name = [:circle,:flower,:parallelogram]
all_params = @strdict weight_approach verbose n ϕ_name order
cases = dict_list(all_params)
β₁ = 20.0
β₂ = 5.0e-1
β₃ = 5.0e-1

# Execute case function
function execute_case(case)
  @unpack weight_approach, verbose, n, ϕ_name, order = case
  case_name = savename(case,"jld2",allowedtypes=(Real, String, Symbol, Function))
  println("Executing case: ",case_name)

  # Case parameters
  if ϕ_name == :circle
    ϕ = level_set(CircleParams(center=center,radius=radius))
  elseif ϕ_name == :flower
    ϕ = level_set(FlowerParams(center=center,radius=radius,n=n_petals))
  elseif ϕ_name == :parallelogram
    ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=1))
  else
    error("Case not recognized")
  end
  if weight_approach == :fraction
    λ = 0.5
  else 
    λ = 1.0
  end
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=(n,n),weight_approach=weight_approach,
    verbose=verbose, order=order,β₁=β₁,β₂=β₂,β₃=β₃,λ=λ)

  # Execute main function
  results = copy(case)
  @timeit to "main_$(case_name)" results["l2norm"], results["cond"] = main_poisson(params)
  results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

  return results
end

# Execute cases
for case in cases
  path = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson")
  filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
  data, file = produce_or_load(path,case,execute_case;filename=filename)
end

# Get data
all_results = collect_results(datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson"))
plot_geom_cases = [:circle,:flower,:parallelogram]
plot_approaches = [:standard,:binary,:fraction]
plot_orders = [1]

# Plot results
labels = ["WSBM","SBM","O-SBM"]
markers = [:circle,:square,:utriangle]
lines = [:solid,:dash]
colors = ["#0072B2", "#E69F00", "#009E73"]
filenames = ["l2norm_circle.pdf","l2norm_flower.pdf","l2norm_parallelogram.pdf"]
shift = 5.0e2
for (igeom, geom_case) in enumerate(plot_geom_cases)
  plt = plot(xlabel="Number of elements",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
  for (iorder, order) in enumerate(plot_orders)
    plot_ref_line = true
    for (iapproach, approach) in enumerate(plot_approaches)
      results = @linq all_results[all_results.:ϕ_name.==geom_case .&& all_results.:weight_approach.==approach .&& all_results.:order.==order,:] |> orderby(:n)
      plot!(plt,results.:n,results.:l2norm,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
      if plot_ref_line
        plot!(plt,results.:n,shift*(results.:n).^(-order-1.),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
        x_triangle = [60, 80, 80]
        y_triangle = [shift*(x_triangle[1]).^(-order-1.),shift*(x_triangle[1]).^(-order-1.),shift*(x_triangle[2]).^(-order-1.)]
        plot!(x_triangle, y_triangle, lw = 1, color = :black, label = "")
        annotate!(70, 1.2*shift*(x_triangle[1]).^(-order-1.), text("1", :black, 9, :left))
        annotate!(83, shift*(70).^(-order-1.), text("$(order+1)", :black, 9, :left))
        plot_ref_line = false
      end
    end
  end
  display(plt)
  savefig(plt,plotsdir("Journal_paper_GSBM","Poisson","ManufacturedPoisson",filenames[igeom]))
end

function plot_vtk_shapes(n)
  n_cells = (n,n)

  # VTK results circle case
  ϕ = level_set(CircleParams(center=center,radius=radius))
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_circle")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

  # VTK results flower case
  ϕ = level_set(FlowerParams(center=center,radius=radius,n=n_petals))
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_flower")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

  # VTK results parallelogram case
  ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=1))
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_parallelogram")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

end


function plot_vtk_approaches(n)
  n_cells = (n,n)
  ϕ = level_set(CircleParams(center=center,radius=radius))

  # VTK results flower WSBM case
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_circle_standard")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,weight_approach=:standard,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

  # VTK results flower SBM case
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_circle_binary")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,weight_approach=:binary,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

  # VTK results flower λ-SBM case
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ManufacturedPoisson","VTK_circle_fraction")
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,weight_approach=:fraction,λ=0.5,β₁=β₁,β₂=β₂,β₃=β₃)
  main_poisson(params)

end

plot_vtk_shapes(40)
plot_vtk_approaches(40)


end

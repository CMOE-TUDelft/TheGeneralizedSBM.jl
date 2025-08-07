module ConditionNumberPoisson
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TheGeneralizedSBM
using Gridap
using TimerOutputs
using Plots
using DataFrames, DataFramesMeta
using Colors

to = TimerOutput("ConditionNumberPoisson")

# Define fixed parameters
center = [0.5,0.5]
radius = 0.3
n_petals = 5
u(x) = sin(15*π*x[1])*sin(15*π*x[2])
f(x) = -tr(∇∇(u)(x))
g(n) = x -> ∇(u)(x)⋅n(x)

# Warm-up parameters
ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=1))
n_cells = (8,8)
output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ConditionNumberPoisson")
params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=n_cells,output_folder=output_folder,order=2)

# Execute main function (Warm-up)
println(main_poisson(params))

# Test parameters
weight_approach = [:fraction]
verbose = false
n = [20,40,80,160]
order = [1]
ϕ_name = [:circle]
radius = collect(0.1:0.1:0.6)
all_params = @strdict weight_approach verbose n ϕ_name order radius
cases = dict_list(all_params)

# Execute case function
function execute_case(case)
  @unpack weight_approach, verbose, n, ϕ_name, order, radius = case
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
  params = PoissonParams(ϕ=ϕ,f=f,g=g,u₀=u,n_cells=(n,n),weight_approach=weight_approach,
    verbose=verbose, order=order,β₁=20.0,β₂=5.0e-1,β₃=5.0e-1,λ=λ,compute_cond=true)

  # Execute main function
  results = copy(case)
  @timeit to "main_$(case_name)" results["l2norm"], results["cond"] = main_poisson(params)
  results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

  return results
end

# Execute cases
for case in cases
  path = datadir("sims","Journal_paper_GSBM","Poisson","ConditionNumberPoisson")
  filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
  data, file = produce_or_load(path,case,execute_case;filename=filename)
end

# Get data
all_results = collect_results(datadir("sims","Journal_paper_GSBM","Poisson","ConditionNumberPoisson"))

# Plot results
# Define the two endpoint colors (e.g., blue and orange for good contrast)
start_color = colorant"#4477AA"  # Blue (colorblind-friendly)
end_color   = colorant"#D55E00"  # Orange (colorblind-friendly)

# Generate 7 interpolated colors
colors = [RGB(
  (1 - t) * start_color.r + t * end_color.r,
  (1 - t) * start_color.g + t * end_color.g,
  (1 - t) * start_color.b + t * end_color.b
) for t in range(0, stop=1, length=6)]

plt = plot(xlabel="Number of elements",ylabel="κ(A)",legend=:topleft,xaxis=:log,yaxis=:log)
plt2 = plot(xlabel="Number of elements",ylabel="Error",legend=:topright,xaxis=:log,yaxis=:log)
global plot_ref_line = true
for (i,radius) in enumerate(radius)
  global plot_ref_line 
  results = @linq all_results[
    all_results.:ϕ_name.==:circle .&& 
    all_results.:weight_approach.==:fraction .&& 
    all_results.:order.==order .&& 
    all_results.:radius.==radius,:] |> orderby(:n)
    plot!(plt,results.:n,results.:cond,marker=:circle,label="R = $radius",color=color=colors[i],xticks = (results.:n, string.(results.:n)),lw=1.5)
    plot!(plt2,results.:n,results.:l2norm,marker=:circle,label="R = $radius",color=color=colors[i],xticks = (results.:n, string.(results.:n)),lw=1.5)
  if plot_ref_line
    shift_4 = 7.0e0
    shift_2 = 5.0e-1
    shift_e = 2.0e2
    plot!(plt,results.:n,shift_4*(results.:n).^(4),xticks = (results.:n, string.(results.:n)),ls=:dash,color=:black,label=false,lw=1.5)
    x_triangle = [40, 40, 50]
    y_triangle = [shift_4*(40).^(4),shift_4*(50).^(4),shift_4*(50).^(4)]
    plot!(x_triangle, y_triangle, lw = 1, color = :black, label = "")
    annotate!(43, 4*shift_4*(40).^(4), text("1", :black, 9, :left))
    annotate!(38, 0.7*shift_4*(50).^(4), text("4", :black, 9, :left))
    plot!(plt,results.:n,shift_2*(results.:n).^(2),ls=:dashdot,color=:black,label=false,lw=1.5)
    x_triangle = [40, 50, 50]
    y_triangle = [shift_2*(40).^(2),shift_2*(40).^(2),shift_2*(50).^(2)]
    plot!(x_triangle, y_triangle, lw = 1, color = :black, label = "")
    annotate!(45, 0.6*shift_2*(40).^(2), text("1", :black, 9, :left))
    annotate!(52, 0.8*shift_2*(50).^(2), text("2", :black, 9, :left))
    plot!(plt2,results.:n,shift_e*(results.:n).^(-2),xticks = (results.:n, string.(results.:n)),ls=:dash,color=:black,label=false,lw=1.5)
    x_triangle = [40, 50, 50]
    y_triangle = [shift_e*(40).^(-2),shift_e*(40).^(-2),shift_e*(50).^(-2)]
    plot!(x_triangle, y_triangle, lw = 1, color = :black, label = "")
    annotate!(45, 1.2*shift_e*(40).^(-2), text("1", :black, 9, :left))
    annotate!(52, 1.2*shift_e*(50).^(-2), text("2", :black, 9, :left))
    global plot_ref_line = false
  end
end        
savefig(plt,plotsdir("Journal_paper_GSBM","Poisson","ConditionNumberPoisson","cond_circle.pdf"))
savefig(plt2,plotsdir("Journal_paper_GSBM","Poisson","ConditionNumberPoisson","l2norm_circle.pdf"))

function plot_solution_vtk(n)
  ϕ = level_set(CircleParams(center=center,radius=0.6))
  output_folder = datadir("sims","Journal_paper_GSBM","Poisson","ConditionNumberPoisson","VTK_radius_0.6")
      params = PoissonParams(ϕ=ϕ,f=f,g=g,u₀=u,n_cells=(n,n),weight_approach=:fraction,
      verbose=true, order=1,β₁=20.0,β₂=5.0e-1,β₃=5.0e-1,λ=0.5,compute_cond=true,output_folder=output_folder)
  main_poisson(params)  
end

plot_solution_vtk(40)
  

end

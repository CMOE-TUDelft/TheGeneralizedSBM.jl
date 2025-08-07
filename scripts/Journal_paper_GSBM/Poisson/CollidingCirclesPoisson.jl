module CollidingCirclesPoisson
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TheGeneralizedSBM
using Gridap
using TimerOutputs
using Plots
using DataFrames, DataFramesMeta

to = TimerOutput("CollidingCirclesPoisson")

# Define fixed parameters
c1 = [0.25,0.5]
c2₀ = [0.75,0.5]
r1 = 0.15
r2 = 0.3
u(x) = sin(15*π*x[1])*sin(15*π*x[2])#s(x[1]+x[2])^2
f(x) = -tr(∇∇(u)(x))

# Warm-up parameters
ϕ1 = level_set(CircleParams(center=c1,radius=r1,in_out=1))
ϕ2 = level_set(CircleParams(center=c2₀,radius=r2,in_out=1))
ϕ =  x -> -1.0*min(ϕ1(x),ϕ2(x))
n_cells = (10,10)
output_folder = datadir("sims","Journal_paper_GSBM","Poisson","CollidingCirclesPoisson")
params = PoissonParams(ϕ=ϕ,n_cells=n_cells,output_folder=output_folder)

# Execute main function (Warm-up)
main_poisson(params)

# Parametric geometry
c1 = [0.2,0.5]
c2(δ) = [0.5+δ,0.5]
r1 = 0.15
r2 = 0.15

# Test parameters
weight_approach = [:standard,:binary,:fraction]
verbose = false
n = [100]#[40,80,160]
δ = collect(-0.2:0.0125:0.2)
all_params = @strdict weight_approach verbose n δ
cases = dict_list(all_params)
β₁ = 20.0
β₂ = 5.0e-1
β₃ = 5.0e-1

# Execute case function
function execute_case(case)
  @unpack weight_approach, verbose, n, δ = case
  case_name = savename(case,"jld2")
  println("Executing case: ",case_name)

  # Case parameters
  ϕ1 = level_set(CircleParams(center=c1,radius=r1,in_out=1))
  ϕ2 = level_set(CircleParams(center=c2(δ),radius=r2,in_out=1))
  ϕ =  x -> -1.0*min(ϕ1(x),ϕ2(x))
  params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=(n,n),weight_approach=weight_approach,
    verbose=verbose,β₁=β₁,β₂=β₂,β₃=β₃)

  # Execute main function
  results = copy(case)
  @timeit to "main_$(case_name)" results["l2norm"], results["cond"] = main_poisson(params)
  results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

  return results
end

# Execute cases
for case in cases
  path = datadir("sims","Journal_paper_GSBM","Poisson","CollidingCirclesPoisson")
  filename = config -> savename(config)
  data, file = produce_or_load(path,case,execute_case;filename=filename)
end

# Get data
n = 100
all_results = collect_results(datadir("sims","Journal_paper_GSBM","Poisson","CollidingCirclesPoisson"))
results_sbm = @linq all_results[all_results.:n.==n.&& all_results.:weight_approach.==:binary,:] |> orderby(:δ)
results_wsbm = @linq all_results[all_results.:n.==n .&& all_results.:weight_approach.==:standard,:] |> orderby(:δ)
l2_sbm = results_sbm.:l2norm
δ_sbm = results_sbm.:δ
t_sbm = results_sbm.:time
l2_wsbm = results_wsbm.:l2norm
δ_wsbm = results_wsbm.:δ
t_wsbm = results_wsbm.:time

# Plot results
plot_approaches = [:standard,:binary,:fraction]
labels = ["WSBM","SBM","λ-SBM"]
markers = [:circle,:square,:utriangle]
lines = [:solid,:dash]
colors = ["#0072B2", "#E69F00", "#009E73"]
plt1 = plot(xlabel="δ",ylabel="L² error",legend=:topleft,yaxis=:log10)
for (iapproach, approach) in enumerate(plot_approaches)
  results = @linq all_results[all_results.:n.==n .&& all_results.:weight_approach.==approach,:] |> orderby(:δ)
  plot!(plt1,results.:δ,results.:l2norm,marker=markers[iapproach],color=colors[iapproach],label=labels[iapproach])
end
# plot!(plt1,δ_wsbm,l2_wsbm,marker=:square,label="WSBM")
plot!(plt1,[0.0], seriestype="vline",ls=:dash,color=:black,label=false)
annotate!(0, 10.0^(-2.45), text(" ⟵ Topology change",Plots.font(10), :left))
# plt3 = plot(xlabel="δ",ylabel="t",legend=:bottomleft)
# plot!(plt3,δ_sbm,t_sbm,marker=:circle,label="SBM")
# plot!(plt3,δ_wsbm,t_wsbm,marker=:square,label="WSBM")
savefig(plt1,plotsdir("Journal_paper_GSBM","Poisson","CollidingCirclesPoisson","l2norm.pdf"))
# savefig(plt3,plotsdir("CollidingCirclesPoisson","time.pdf"))

# VTK results cases
function plot_vtk(n)
  weight_approach = :standard
  verbose = true
  δs = [0.15,0.175,0.2,-0.1,0.0,0.1]
  for δ in δs
    ϕ1 = level_set(CircleParams(center=c1,radius=r1,in_out=1))
    ϕ2 = level_set(CircleParams(center=c2(δ),radius=r2,in_out=1))
    ϕ =  x -> -1.0*min(ϕ1(x),ϕ2(x))
    output_folder = datadir("sims","Journal_paper_GSBM","Poisson","CollidingCirclesPoisson","VTK_delta_$δ")
    params = PoissonParams(ϕ=ϕ,f=f,u₀=u,n_cells=(n,n),weight_approach=weight_approach,
      verbose=verbose,output_folder=output_folder,β₁=β₁,β₂=β₂,β₃=β₃)
      main_poisson(params)
  end
end

plot_vtk(n)

end

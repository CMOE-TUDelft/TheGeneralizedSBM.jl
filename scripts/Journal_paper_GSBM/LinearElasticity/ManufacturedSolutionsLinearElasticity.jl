function run_manufactured_solutions_linear_elasticity()

  to = TimerOutput("ManufacturedLinearElasticity")

  # Define fixed parameters
  center = [0.5,0.5]
  radius = 0.3
  n_petals = 5
  a₀ = 1/10
  u₁(x) = -a₀*cos(3π*x[1])*sin(π*x[2])
  u₂(x) = a₀*sin(7π*x[1])*sin(5π*x[2])
  u(x) = VectorValue(u₁(x),u₂(x))
  E = 200.0e9
  ν = 0.3
  λₑ = ν*E/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  I = TensorValue(1.0,0.0,0.0,1.0)
  ∂u₁∂x₁(x) = a₀*3π*sin(3π*x[1])*sin(π*x[2])
  ∂u₁∂x₂(x) = -a₀*π*cos(3π*x[1])*cos(π*x[2])
  ∂u₂∂x₁(x) = a₀*7π*cos(7π*x[1])*sin(5π*x[2])
  ∂u₂∂x₂(x) = a₀*5π*sin(7π*x[1])*cos(5π*x[2])
  f(x) = -(λₑ+μ)*(∇(∂u₁∂x₁)(x)+∇(∂u₂∂x₂)(x)) - μ*VectorValue((∇(∂u₁∂x₁)(x))⋅VectorValue(1,0)+(∇(∂u₁∂x₂)(x))⋅VectorValue(0,1),
  (∇(∂u₂∂x₁)(x))⋅VectorValue(1,0)+(∇(∂u₂∂x₂)(x))⋅VectorValue(0,1))

  # Warm-up parameters
  ϕ = level_set(CircleParams(center=center,radius=radius))
  n_cells = (10,10)
  output_folder = datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity")
  params = LinearElasticityParams(ϕ=ϕ,f=f,u₀=u,λₑ=λₑ,μ=μ,n_cells=n_cells,output_folder=output_folder)

  # Execute main function (Warm-up)
  println(main_elasticity(params))

  # Test parameters
  weight_approach = [:standard,:binary,:fraction]
  verbose = false
  n = [20,40,80,160,320]
  order = [1]#,2]
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
    params = LinearElasticityParams(ϕ=ϕ,f=f,u₀=u,λₑ=λₑ,μ=μ,n_cells=(n,n),weight_approach=weight_approach,
      verbose=verbose,order=order,β₁=β₁,β₂=β₂,β₃=β₃,λ=λ)

    # Execute main function
    results = copy(case)
    @timeit to "main_$(case_name)" results["l2norm"] = main_elasticity(params)
    results["time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

    return results
  end

  # Execute cases
  for case in cases
    path = datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity")
    filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
    data, file = produce_or_load(path,case,execute_case;filename=filename)
  end

  # Get data
  all_results = collect_results(datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity"))
  plot_geom_cases = [:circle,:flower,:parallelogram]
  plot_approaches = [:standard,:binary,:fraction]
  plot_orders = [1]#,2]

  # Plot results
  labels = ["WSBM","SBM","λ-SBM"]
  markers = [:circle,:square,:utriangle]
  lines = [:solid,:dash]
  colors = ["#0072B2", "#E69F00", "#009E73"]
  filenames = ["l2norm_circle.pdf","l2norm_flower.pdf","l2norm_parallelogram.pdf"]
  shift = 1.0e1
  for (igeom, geom_case) in enumerate(plot_geom_cases)
    plt = plot(xlabel="Number of elements",ylabel="L² error",legend=:topright,xaxis=:log,yaxis=:log)
    plot_ref_line = true
    for (iorder, order) in enumerate(plot_orders)
      for (iapproach, approach) in enumerate(plot_approaches)
        results = @linq all_results[all_results.:ϕ_name.==geom_case .&& all_results.:weight_approach.==approach .&& all_results.:order.==order,:] |> orderby(:n)
        plot!(plt,results.:n,results.:l2norm,marker=markers[iapproach],ls=lines[iorder],color=colors[iapproach],label=labels[iapproach])
        if plot_ref_line
          # plot!(plt,1.0./results.:n,2*(1.0./results.:n).^(1),ls=:dash,color=:black,label=false)
          plot!(plt,results.:n,shift*(results.:n).^(-2),xticks = (results.:n, string.(results.:n)),ls=:dashdot,color=:black,label=false)
          # plot!(plt,1.0./results.:n,1.0e-1*(1.0./results.:n).^(3),ls=:dot,color=:black,label=false)
          x_triangle = [80, 100, 100]
          y_triangle = [shift*(80).^(-2),shift*(80).^(-2),shift*(100).^(-2)]
          plot!(x_triangle, y_triangle, lw = 1, color = :black, label = "")
          annotate!(90, 1.2*shift*(80).^(-2), text("1", :black, 9, :left))
          annotate!(103, shift*(90).^(-2), text("2", :black, 9, :left))
          plot_ref_line = false
        end
      end
    end
    savefig(plt,plotsdir("Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity",filenames[igeom]))
  end

  function plot_vtk_shapes(n)
    n_cells = (n,n)

    # VTK results circle case
    ϕ = level_set(CircleParams(center=center,radius=radius))
    output_folder = datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity","VTK_circle")
    params = LinearElasticityParams(ϕ=ϕ,f=f,u₀=u,λₑ=λₑ,μ=μ,n_cells=n_cells,
      output_folder=output_folder,weight_approach=:fraction,λ=0.5,β₁=β₁,β₂=β₂,β₃=β₃)
    main_elasticity(params)

    # VTK results flower case
    ϕ = level_set(FlowerParams(center=center,radius=radius,n=n_petals))
    output_folder = datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity","VTK_flower")
    params = LinearElasticityParams(ϕ=ϕ,f=f,u₀=u,λₑ=λₑ,μ=μ,n_cells=n_cells,
      output_folder=output_folder,weight_approach=:fraction,λ=0.5,β₁=β₁,β₂=β₂,β₃=β₃)
    main_elasticity(params)

    # VTK results parallelogram case
    ϕ = level_set(ParallelogramParams(v₁=[0.21,0.21],v₂=[0.64,0.31],v₃=[0.83,0.83],v₄=[0.37,0.73],in_out=1))
    output_folder = datadir("sims","Journal_paper_GSBM","LinearElasticity","ManufacturedLinearElasticity","VTK_parallelogram")
    params = LinearElasticityParams(ϕ=ϕ,f=f,u₀=u,λₑ=λₑ,μ=μ,n_cells=n_cells,
      output_folder=output_folder,weight_approach=:fraction,λ=0.5,β₁=β₁,β₂=β₂,β₃=β₃)
    main_elasticity(params)
  end

  plot_vtk_shapes(40)

end

function run_cylinder_cross_flow()

  to = TimerOutput("CylinderCrossFlow")

  # Define geometry parameters
  domain = (-10.0,5.0,-6.0,6.0)
  y₀ = 0.5
  f₀ = 0.2237
  y(t) = y₀*cos(2π*f₀*t)
  center(t) = [-5.0,y(t)]
  radius = 1.0

  # Define the velocity field and forcing terms
  uy₀ = 2π*f₀*y₀
  ucy(t) = -uy₀*sin(2π*f₀*t)
  uc(x,t) = VectorValue(ucy(t),0.0)
  uc(t::Real) = x -> uc(x,t)
  uin(t::Real) = x -> VectorValue(1.0,0.0)
  f(t::Real) = x -> VectorValue(0.0,0.0)
  g(t::Real) = x -> 0.0

  # Define algorithmic parameters
  β₁=1.0
  βᵤ=1.0
  βₚ=1.0e-1
  βdiv=1.0
  order = 2
  ν = 0.01

  # Define time step size
  h = 0.2
  @show ΔtCFL_D_1 = 2*radius/uy₀
  @show ΔtCFL_h_1 = h/uy₀

  # Warm-up parameters
  ϕ(t) = level_set(CircleParams(center=center(t),radius=radius,in_out=1))
  n_cells = (11,7)
  output_folder = datadir("sims","Journal_paper_GSBM","Stokes","CylinderCrossFlow")
  ode_solver_params = TimeIntegratorParams(T=0.1,method=:generalized_alpha,ρ∞=0.0)
  params = TransientStokesParams(
    domain=domain,
    ϕ=ϕ,
    u₀=uc,
    ν=ν,
    n_cells=n_cells,
    output_folder=output_folder,
    order=order,
    ode_solver_params=ode_solver_params,
    is_exact_solution=false,
    in_out_wall_tags=([1,3,7],[8],[2,4,5,6]),
    dirichlet_tags=["inlet","wall"],
    dirichlet_masks=[(true,true),(false,true)],
    dirichlet_functions=[uin,uin],
    verbose=true,
  )

  # Execute main function (Warm-up)
  # println(main_transient_stokes(params))

  # CFL test parameters
  weight_approach = [:standard]
  CFL = [1.0,0.5,0.1]
  global_gp = [true,false]
  all_params = @strdict weight_approach CFL global_gp
  cases = dict_list(all_params)

  # Execute case function
  function execute_case(case)
    @unpack weight_approach, CFL, global_gp = case
    case_name = savename(case,"jld2",allowedtypes=(Real, String, Symbol, Function))
    println("Executing case: $case_name")

    # Define parameters
    L = domain[2] - domain[1]
    @show n = ceil(Int,L/h)
    @show Δt = CFL*h/uy₀
    @show T = 2/f₀
    ode_solver_params = TimeIntegratorParams(Δt=Δt,T=T,method=:generalized_alpha,ρ∞=0.0)
    verbose = true
    output_folder = datadir("sims","Journal_paper_GSBM","Stokes","CylinderCrossFlow",replace(case_name,".jld2"=>""))
    if !isdir(output_folder)
      mkpath(output_folder)
    end

    # Case parameters
    params = TransientStokesParams(
      domain=domain,
      ϕ=ϕ,
      u₀=uc,
      ν=ν,
      n_cells=(n,ceil(Int,n/2)),
      output_folder=output_folder,
      order=order,
      weight_approach=weight_approach,
      global_gp=global_gp,
      ode_solver_params=ode_solver_params,
      is_exact_solution=false,
      in_out_wall_tags=([1,3,7],[8],[2,4,5,6]),
      dirichlet_tags=["inlet","wall"],
      dirichlet_masks=[(true,true),(false,true)],
      dirichlet_functions=[uin,uin],
      β₁=β₁,
      βᵤ=βᵤ,
      βₚ=βₚ,
      βdiv=βdiv,
      verbose=true,
    )
    println("Parameters: ",params)

    # Execute main function
    results = copy(case)
    @timeit to "main_$(case_name)" results["l2l2ᵤ"], results["l2h1ᵤ"], results["l2l2ₚ"], results["l2h1ₚ"], results["l∞l2ᵤ"], results["l∞h1ᵤ"], results["l∞l2ₚ"], results["l∞h1ₚ"], results["lastl2ᵤ"], results["lastl2ₚ"], results["FD"], results["time"] = main_transient_stokes(params)
    results["cpu_time"] = TimerOutputs.time(to["main_$(case_name)"])/1.0e9

    return results

  end

  # Execute CFL cases
  for case in cases
    path = datadir("sims","Journal_paper_GSBM","Stokes","CylinderCrossFlow")
    filename = config -> savename(config,allowedtypes=(Real, String, Symbol, Function))
    data, file = produce_or_load(path,case,execute_case;filename=filename)
  end

  # Get data
  all_results = collect_results(datadir("sims","Journal_paper_GSBM","Stokes","CylinderCrossFlow"))

  
  colors = ["#0072B2", "#E69F00", "#009E73"]
  xlims = (1.0,8.0)
  for iglobal_gp in global_gp
    println("Plotting results for global_gp=$(iglobal_gp)")
    plt_name = "global_gp=$(iglobal_gp)"
    plt = plot(xlabel="Time",ylabel="Pressure",legend=:bottomright,xlims=xlims,lw=2)
    plt2 = plot(xlabel="Time",ylabel="FD",legend=:topright,xlims=xlims)
    for (i,iCFL) in enumerate(CFL)
      println("  CFL=$(iCFL)")
      data = CSV.File(datadir("sims","Journal_paper_GSBM","Stokes","CylinderCrossFlow","pressure_probe_CFL=$(iCFL)_global_gp=$(iglobal_gp).csv"))
      results = @linq all_results[all_results.:global_gp.==iglobal_gp .&& all_results.:CFL.==iCFL ,:]
      # println("CFL=$(iCFL), global_gp=$(iglobal_gp): ",results)
      t = data["Time"]
      p = data["pₕ"]
      plot!(plt,t,p,ls=:solid,color=colors[i],label="CFL=$(iCFL)",lw=1.2)
      plot!(plt2,results.:time,results.:FD,ls=:solid,color=colors[i],label="CFL=$(iCFL)",lw=1.2)
    end
    savefig(plt,plotsdir("Journal_paper_GSBM","Stokes","CylinderCrossFlow",plt_name*".pdf"))
    savefig(plt2,plotsdir("Journal_paper_GSBM","Stokes","CylinderCrossFlow",plt_name*"_FD.pdf"))
  end
  return nothing
end
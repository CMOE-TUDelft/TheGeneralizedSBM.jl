module TheGeneralizedSBM
using DrWatson
@quickactivate "TheGeneralizedSBM"
using TimerOutputs
using Plots
using Colors
using DataFrames, DataFramesMeta
using Parameters
using LinearAlgebra
using Gridap
using GridapSolvers
using GridapSolvers.LinearSolvers

include("HelperFunctions.jl")
include("TimeIntegrator.jl")
include("Poisson.jl")
include("LinearElasticity.jl")
include("Stokes.jl")
include("TransientStokes.jl")

export compute_weights
export level_set
export CircleParams, FlowerParams, ParallelogramParams, SphereParams
export PoissonParams, main_poisson
export LinearElasticityParams, main_elasticity
export StokesParams, main_stokes
export NavierStokesParams, main_navier_stokes
export TransientStokesParams, main_transient_stokes
export TimeIntegratorParams

# Include script files
include("../scripts/Journal_paper_GSBM/Poisson/ManufacturedSolutionsPoisson.jl")
include("../scripts/Journal_paper_GSBM/Poisson/ConditionNumberPoisson.jl")
include("../scripts/Journal_paper_GSBM/Poisson/CollidingCirclesPoisson.jl")
include("../scripts/Journal_paper_GSBM/LinearElasticity/ManufacturedSolutionsLinearElasticity.jl")
include("../scripts/Journal_paper_GSBM/Stokes/ManufacturedSolutionsStokes.jl")
include("../scripts/Journal_paper_GSBM/Stokes/ManufacturedSolutionsTransientStokes_h.jl")
include("../scripts/Journal_paper_GSBM/Stokes/ManufacturedSolutionsTransientStokes_dt.jl")
include("../scripts/Journal_paper_GSBM/Stokes/ManufacturedSolutionsTransientStokes_betap.jl")

function run_tests(test_type::String)
  if test_type == "all"
    # Poisson tests
    run_manufactured_solutions_poisson()
    run_condition_number_poisson()
    run_colliding_circles_poisson()
    # Linear Elasticity tests
    run_manufactured_solutions_linear_elasticity()
    # Stokes tests
    run_manufactured_solutions_stokes()
    run_manufactured_solutions_transient_stokes_h()
    run_manufactured_solutions_transient_stokes_dt()
    run_manufactured_solutions_transient_stokes_betap()
  elseif test_type == "ManufacturedSolutionsPoisson.jl"
    run_manufactured_solutions_poisson()
  elseif test_type == "ConditionNumberPoisson.jl"
    run_condition_number_poisson()
  elseif test_type == "CollidingCirclesPoisson.jl"
    run_colliding_circles_poisson()
  elseif test_type == "ManufacturedSolutionsLinearElasticity.jl"
    run_manufactured_solutions_linear_elasticity()
  elseif test_type == "ManufacturedSolutionsStokes.jl"
    run_manufactured_solutions_stokes()
  elseif test_type == "ManufacturedSolutionsTransientStokes_h.jl"
    run_manufactured_solutions_transient_stokes_h()
  elseif test_type == "ManufacturedSolutionsTransientStokes_dt.jl"
    run_manufactured_solutions_transient_stokes_dt()
  elseif test_type == "ManufacturedSolutionsTransientStokes_betap.jl"
    run_manufactured_solutions_transient_stokes_betap()
  else
    error("Unknown test type: $test_type")
  end
end

end # module TheGeneralizedSBM

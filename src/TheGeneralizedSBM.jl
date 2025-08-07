module TheGeneralizedSBM
using Gridap
using Parameters
using DrWatson
using LinearAlgebra
using GridapSolvers
using GridapSolvers.LinearSolvers

include("HelperFunctions.jl")
include("TimeIntegrator.jl")
include("Poisson.jl")
include("LinearElasticity.jl")
include("Stokes.jl")
include("NavierStokes.jl")
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

end # module TheGeneralizedSBM

"""
    TimeIntegratorParams

Parameters for time integrators.
"""
@with_kw struct TimeIntegratorParams
  Δt::Float64 = 0.1 # Time step
  T::Float64 = 1.0 # Final time
  method::Symbol = :theta_method # Time integration method  
  θ::Float64 = 0.5 # θ parameter for theta method
  ρ∞::Float64 = 1.0 # ρ∞ parameter for generalized alpha method
end

function get_ode_solver(params::TimeIntegratorParams,ls)
  @unpack Δt,T,method,θ,ρ∞ = params
  if method == :theta_method
    return ThetaMethod(ls,Δt,θ)
  elseif method == :generalized_alpha
    return GeneralizedAlpha1(ls,Δt,ρ∞)
  else
    error("ODE solver method not implemented")
  end
end

function get_initial_condition(params::TimeIntegratorParams,x₀,v₀)
  @unpack Δt,T,method,θ = params
  if method == :theta_method
    return x₀
  elseif method == :generalized_alpha
    return (x₀,v₀)
  else
    error("ODE solver method not implemented")
  end
end
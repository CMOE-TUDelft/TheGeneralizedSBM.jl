@with_kw struct CircleParams
  case::Symbol = :circle
  center::Vector{Float64} = [0.0,0.0]
  radius::Float64 = 1.0
  in_out::Int = -1 # -1 if inside, 1 if outside
end

@with_kw struct FlowerParams
  case::Symbol = :flower
  center::Vector{Float64} = [0.0,0.0]
  radius::Float64 = 1.0
  in_out::Int = -1 # -1 if inside, 1 if outside
  a::Float64 = 0.3 # Amplitude of flower oscillations
  n::Int = 3 # number of flower petals
end

@with_kw struct ParallelogramParams
  case::Symbol = :parallelogram
  v₁::Vector{Float64} = [0.0,0.0]
  v₂::Vector{Float64} = [1.0,0.0]
  v₃::Vector{Float64} = [1.0,1.0]
  v₄::Vector{Float64} = [0.0,1.0]
  in_out::Int = 1 # 1 if inside, -1 if outside
end

@with_kw struct SphereParams
  case::Symbol = :sphere
  center::Vector{Float64} = [0.0,0.0,0.0]
  radius::Float64 = 1.0
  in_out::Int = -1 # -1 if inside, 1 if outside
end

"""
    level_set(params)

Compute the level set function `ϕ` from a set of parameters `params`. Currently,
the following cases are implemented:

- `:circle`: a circle with center `center` and radius `radius` and `in_out` = -1 if
  inside and 1 if outside.
"""
function level_set(params)
  @unpack case = params
  if case == :circle
    @unpack center,radius,in_out = params
    return x -> in_out*(√((x[1]-center[1])^2 + (x[2]-center[2])^2) - radius)
  elseif case == :flower
    @unpack center,radius,a,n,in_out = params
    w(x) = VectorValue(x[1]-center[1],x[2]-center[2])
    t(x) = angle(Complex(w(x)[1],w(x)[2]))
    R(x) = radius*(1.0 + a*sin(n*t(x)))
    return x -> in_out*(√(w(x)⋅w(x)) - R(x))
  elseif case == :parallelogram
    @unpack v₁,v₂,v₃,v₄,in_out = params
    x₁,y₁ = v₁
    x₂,y₂ = v₂
    x₃,y₃ = v₃
    x₄,y₄ = v₄
    l₁(x) = (x[2]-y₁)*(x₂ - x₁) - (x[1]-x₁)*(y₂ - y₁)
    l₂(x) = (x[2]-y₂)*(x₃ - x₂) - (x[1]-x₂)*(y₃ - y₂)
    l₃(x) = (x[2]-y₃)*(x₄ - x₃) - (x[1]-x₃)*(y₄ - y₃)
    l₄(x) = (x[2]-y₄)*(x₁ - x₄) - (x[1]-x₄)*(y₁ - y₄)
    d₁(x) = l₁(x)/√((x₂ - x₁)^2 + (y₂ - y₁)^2)
    d₂(x) = l₂(x)/√((x₃ - x₂)^2 + (y₃ - y₂)^2)
    d₃(x) = l₃(x)/√((x₄ - x₃)^2 + (y₄ - y₃)^2)
    d₄(x) = l₄(x)/√((x₁ - x₄)^2 + (y₁ - y₄)^2)
    return x -> in_out*minimum([d₁(x),d₂(x),d₃(x),d₄(x)])
  elseif case == :sphere
    @unpack center,radius,in_out = params
    return x -> in_out*(√((x[1]-center[1])^2 + (x[2]-center[2])^2 + (x[3]-center[3])^2) - radius)
  else
    return error("Case not recognized")
  end
end

"""
    compute_weights(Ω,params)

Compute the weights `α` from a level set function `ϕ`. The weights are computed
using either a binary approach, i.e., α = 1 if ∫ₖ(ϕ)dK /∫ₖ(1)dK ==1 and α = 0 otherwise,
or a standard approach, i.e., α = max(∫ₖ(ϕ)dK /∫ₖ(1)dK,weight_min_value)
"""
function compute_weights(Ω::Triangulation,params)
  @unpack ϕ, weight_quad_degree, weight_min_value, weight_approach = params
  quad = CellQuadrature(Ω,weight_quad_degree)
  if weight_approach == :binary
    return α = CellField(float((integrate(x->(ϕ(x)>0.0),quad)./integrate(1.0,quad)) .≈ 1.0),Ω)
  elseif weight_approach == :fraction
    @unpack λ = params
    return α = CellField(float((integrate(x->(ϕ(x)>0.0),quad)./integrate(1.0,quad)) .>= (1.0-λ)),Ω)
  elseif weight_approach == :standard
     α = CellField(max.(integrate(x->(ϕ(x)>0.0),quad)./integrate(1.0,quad),weight_min_value),Ω)
  else
    return error("Approach not recognized")
  end
end

"""
    d(params)

Compute the distance field from a level set function `ϕ`.
"""
function d(params)
  @unpack ϕ = params
  norm_∇ϕ(x) = max(norm(∇(ϕ)(x)),1.0e-18)
  return x->(-(ϕ(x)))*∇(ϕ)(x)/norm_∇ϕ(x)
end

"""
    meanᵧ(α,v)

Compute the mean of a field `v` using the weights `α`. The mean is computed as
γ*v.⁺ + (1-γ)*v.⁻, with γ = α.⁺/(α.⁺+α.⁻).
"""
function meanᵧ(α,v)
  γ = α.⁺/max.(α.⁺+α.⁻,1.0e-14)
  return γ*v.⁺ + (1-γ)*v.⁻
end

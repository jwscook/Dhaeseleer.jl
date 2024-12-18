module Dhaeseleer

using Symbolics, LinearAlgebra

sqrtrule1 = @rule sqrt((~x)^2) => ~x
sqrtrule2 = @rule sqrt(1/(~x)^2) => 1/(~x)

function process(a)
  a = Symbolics.expand_derivatives(a, true)
  a = Symbolics.expand(a)
  a = Symbolics.simplify(a, Rewriters.Prewalk(Rewriters.PassThrough(sqrtrule1)))
  a = Symbolics.simplify(a, Rewriters.Prewalk(Rewriters.PassThrough(sqrtrule2)))
  #a = Symbolics.simplify(a)
  return a
end

abstract type AbstractCoordinateSystem end

struct ArbitraryCoordinateSystem <: AbstractCoordinateSystem
  xi::Vector{Num}
  ui::Vector{Num}
  ∂₁::Differential
  ∂₂::Differential
  ∂₃::Differential
  ∂u⃗_∂x⃗::Matrix{Num} # ∂ui/∂xj
  gᵢⱼ::Matrix{Num}
  gⁱʲ::Matrix{Num}
  J::Num
end

#gᵢⱼ(c::AbstractCoordinateSystem) = c.gᵢⱼ # c.∂u⃗_∂x⃗ * c.∂u⃗_∂x⃗' #
#gⁱʲ(c::AbstractCoordinateSystem) = c.gⁱʲ # inv(gᵢⱼ(c))        #
#jac(c::AbstractCoordinateSystem) = c.J   # sqrt(det(gᵢⱼ(c)))  #
gⁱʲ(c::AbstractCoordinateSystem) = c.∂u⃗_∂x⃗ * c.∂u⃗_∂x⃗' # c.gⁱʲ #
gᵢⱼ(c::AbstractCoordinateSystem) = inv(gⁱʲ(c))        # c.gᵢⱼ #
jac(c::AbstractCoordinateSystem) = sqrt(det(gᵢⱼ(c)))  # c.J   #

abstract type AbstractDifferentialOperator end
for f in (:gⁱʲ, :gᵢⱼ, :jac)
  @eval $f(d::AbstractDifferentialOperator) = $f(coordinatesystem(d))
end

struct ∇ <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end
struct ∇o <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end
struct ∇x <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end
coordinatesystem(a::AbstractDifferentialOperator) = a.c

function CoordinateSystem(xi, ui)
  ∂₁ = Differential(ui[1])
  ∂₂ = Differential(ui[2])
  ∂₃ = Differential(ui[3])
  ∂u⃗_∂x⃗ = Matrix{Num}(undef, 3, 3)
  ss(x) = String(Symbol(x))
  #for (i, u) in enumerate(ui), (j, x) in enumerate(xi)
  #  ∂u⃗_∂x⃗[i, j] = Symbolics.variable("∂" * ss(u) * "_∂" * ss(x))
  #end
  ∂x = Differential(xi[1])
  ∂y = Differential(xi[2])
  ∂z = Differential(xi[3])
  ∂u⃗_∂x⃗ = [∂x(ui[1]) ∂y(ui[1]) ∂z(ui[1]);
           ∂x(ui[2]) ∂y(ui[2]) ∂z(ui[2]);
           ∂x(ui[3]) ∂y(ui[3]) ∂z(ui[3])]
  #∂u⃗_∂x⃗ .= process.(∂u⃗_∂x⃗)
  @variables g¹¹ g¹² g¹³ g²¹ g²² g²³ g³¹ g³² g³³
  gⁱʲ = [g¹¹ g¹² g¹³; g²¹ g²² g²³; g³¹ g³² g³³] # gⁱʲ = inv(gᵢⱼ)
  @variables g₁₁ g₁₂ g₁₃ g₂₁ g₂₂ g₂₃ g₃₁ g₃₂ g₃₃
  gᵢⱼ = [g₁₁ g₁₂ g₁₃; g₂₁ g₂₂ g₂₃; g₃₁ g₃₂ g₃₃] # ∂u⃗_∂x⃗ * ∂u⃗_∂x⃗'
  @variables J # sqrt(det(gᵢⱼ))
  c = ArbitraryCoordinateSystem(xi, ui, ∂₁, ∂₂, ∂₃, ∂u⃗_∂x⃗, gᵢⱼ, gⁱʲ, J)
  return (∇(c), ∇o(c), ∇x(c))
end

function ArbitraryCoordinateSystem()
  @variables x y z # cartesian
  @variables u1(x, y, z) u2(x, y, z) u3(x, y, z)
  return ArbitraryCoordinateSystem([x, y, z], [u1, u2, u3])
end

abstract type AbstractCoconVector <: AbstractArray{Num, 1} end
for (CV, Ai) in ((:CovariantVector, :Aᵢ), (:ContravariantVector, :Aⁱ))
  @eval struct $(CV){C<:AbstractCoordinateSystem} <: AbstractCoconVector
    cs::C
    $(Ai)::Vector{Num} # for a component Aⁱ uᵢ, this represents Aⁱ
    function $(CV)(cs::C, Ai::Vector{Num}) where {C<:AbstractCoordinateSystem}
      return new{C}(cs, Ai)
    end
  end
  @eval Base.setindex!(a::$(CV), v, i) = (a.$(Ai)[i] = v)
  @eval Base.getindex(a::$(CV), i) = a.$(Ai)[i]
  @eval function Base.iterate(a::$(CV), state=1)
    return state <= 3 ? (a.$(Ai)[state], state+1) : nothing
  end
end
Base.eachindex(::AbstractCoconVector) = 1:3
Base.length(::AbstractCoconVector) = 3
Base.size(::AbstractCoconVector) = (3,)
function simplify!(input::AbstractArray)
  for i in eachindex(input)
    input[i] = process(input[i])
  end
  return input
end
simplify(input::Num) = process(input)
function subsimp!(input::AbstractArray, dict)
  for i in eachindex(input)
    input[i] = Symbolics.substitute(input[i], dict)
  end
  return simplify!(input)
end
subsimp(input::Num, dict) = simplify(Symbolics.substitute(input, dict))

coordinatesystem(a) = a.cs

function Base.convert(a::CovariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, gⁱʲ(cs) * a.Aᵢ)
end

function Base.convert(a::ContravariantVector)
  cs = coordinatesystem(a)
  return CovariantVector(cs, gᵢⱼ(cs) * a.Aⁱ)
end

function normalise(a::ContravariantVector)
  cs = coordinatesystem(a)
  #Aⁱ uᵢ = Ai ei
  #ei = uᵢ / |uᵢ|
  #Aⁱ uᵢ = Ai uᵢ / |uᵢ|
  #Aⁱ |uᵢ| = Ai
  gij = gᵢⱼ(cs)
  return [a.Aⁱ[1] * sqrt(gij[1, 1]),
          a.Aⁱ[2] * sqrt(gij[2, 2]),
          a.Aⁱ[3] * sqrt(gij[3, 3])]
end

function normalise(a::CovariantVector)
  cs = coordinatesystem(a)
  #Aᵢ uⁱ = Ai ei
  #ei = uⁱ / |uⁱ|
  #Aᵢ uⁱ = Ai uⁱ / |uⁱ|
  #Aᵢ |uⁱ| = Ai
  gij = gⁱʲ(cs)
  return [a.Aᵢ[1] * sqrt(gij[1, 1]),
          a.Aᵢ[2] * sqrt(gij[2, 2]),
          a.Aᵢ[3] * sqrt(gij[3, 3])]
end

import Base: abs2, abs
# aᵢ uⁱ = a_ψ ∇ψ + a_θ ∇θ + a_φ ∇φ
abs(a::AbstractCoconVector) = sqrt(abs2(a))
Base.:*(a::AbstractCoconVector, b::Num) = b * a
for (CV, Ai) in ((CovariantVector, :Aᵢ), (ContravariantVector, :Aⁱ))
  for op in (:+, :-)
    @eval function Base.$(op)(a::T, b::T) where {T<:$(CV)}
      cs = coordinatesystem(a)
      @assert cs == coordinatesystem(b)
      return CovariantVector(cs, $(op)(a.$(Ai), b.$(Ai)))
    end
  end
  for op in (:*, :/)
    @eval function Base.$(op)(a::$(CV), b::Num)
      return CovariantVector(coordinatesystem(a), $(op).(a.$(Ai), b))
    end
  end
  @eval Base.:-(a::$(CV)) = CovariantVector(coordinatesystem(a), -a.$(Ai))
  @eval abs2(a::$(CV)) = (a.$(Ai)' * a.$(Ai))
end
Base.:*(a::Num, b::AbstractCoconVector) = b * a
Base.:+(a::CovariantVector, b::ContravariantVector) = convert(a) + b
Base.:+(a::ContravariantVector, b::CovariantVector) = convert(a) + b
Base.:-(a::CovariantVector, b::ContravariantVector) = convert(a) - b
Base.:-(a::ContravariantVector, b::CovariantVector) = convert(a) - b

function (d::∇)(q)
  return CovariantVector(d.c, [d.c.∂₁(q), d.c.∂₂(q), d.c.∂₃(q)])
end

function (d::∇o)(a::ContravariantVector)
  cs = coordinatesystem(a)
  J = jac(cs)
  f = (cs.∂₁(J * a.Aⁱ[1]) + cs.∂₂(J * a.Aⁱ[2]) + cs.∂₃(J * a.Aⁱ[3])) / J
  return f
end
∇o(a::CovariantVector) = ∇o(convert(a))

function (d::∇x)(a::CovariantVector)
  cs = coordinatesystem(a)
  J = jac(cs)
  return ContravariantVector(cs, [(cs.∂₂(a.Aᵢ[3]) - cs.∂₃(a.Aᵢ[2])) / J,
                                  (cs.∂₃(a.Aᵢ[1]) - cs.∂₁(a.Aᵢ[3])) / J,
                                  (cs.∂₁(a.Aᵢ[2]) - cs.∂₂(a.Aᵢ[1])) / J])
end
(d::∇x)(a::ContravariantVector) = d(convert(a))

import LinearAlgebra: cross, dot
dot(a::ContravariantVector, b::CovariantVector) = dot(a.Aⁱ, b.Aᵢ)
dot(a::CovariantVector, b::ContravariantVector) = dot(b, a)
dot(a::T, b::T) where {T<:AbstractCoconVector} = dot(a, convert(b))
function cross(a::CovariantVector, b::CovariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, cross(a.Aᵢ, b.Aᵢ) ./ jac(cs))
end
function cross(a::ContravariantVector, b::ContravariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, cross(a.Aⁱ, b.Aⁱ) .* jac(cs))
end
cross(a::CovariantVector, b::ContravariantVector) = cross(convert(a), b)
cross(a::ContravariantVector, b::CovariantVector) = cross(a, convert(b))
export cross, dot

end # module


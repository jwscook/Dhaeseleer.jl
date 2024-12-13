module Dhaeseleer

using Symbolics, LinearAlgebra

abstract type AbstractCoordinateSystem end

struct ArbitraryCoordinateSystem <: AbstractCoordinateSystem
  xi::Vector{Num}
  ui::Vector{Num}
  ∂₁::Differential
  ∂₂::Differential
  ∂₃::Differential
  ∂u⃗_∂x⃗::Matrix{Num} # ∂ui/∂xj
end

struct ∇
  c::AbstractCoordinateSystem
end
struct ∇o
  c::AbstractCoordinateSystem
end
struct ∇x
  c::AbstractCoordinateSystem
end

function CoordinateSystem(xi, ui)
  ∂x = Differential(xi[1])
  ∂y = Differential(xi[2])
  ∂z = Differential(xi[3])
  ∂₁ = Differential(ui[1])
  ∂₂ = Differential(ui[2])
  ∂₃ = Differential(ui[3])
  ∂u⃗_∂x⃗ = [∂x(ui[1]) ∂y(ui[1]) ∂z(ui[1]);
           ∂x(ui[2]) ∂y(ui[2]) ∂z(ui[2]);
           ∂x(ui[3]) ∂y(ui[3]) ∂z(ui[3])]
  c = ArbitraryCoordinateSystem(xi, ui, ∂₁, ∂₂, ∂₃, ∂u⃗_∂x⃗)
  return (∇(c), ∇o(c), ∇x(c))
end

function ArbitraryCoordinateSystem()
  @variables x y z # cartesian
  @variables u1(x, y, z) u2(x, y, z) u3(x, y, z)
  return ArbitraryCoordinateSystem([x, y, z], [u1, u2, u3])
end

abstract type AbstractCoconVector end
for (CV, Ai) in ((:CovariantVector, :Aᵢ), (:ContravariantVector, :Aⁱ))
  @eval struct $(CV){C<:AbstractCoordinateSystem} <: AbstractCoconVector
    cs::C
    $(Ai)::Vector{Num} # for a component Aⁱ uᵢ, this represents Aⁱ
    function $(CV)(cs::C, Ai::Vector{Num}) where {C<:AbstractCoordinateSystem}
      return new{C}(cs, Symbolics.expand_derivatives(Ai))
    end
  end
end

#struct CovariantVector{C<:AbstractCoordinateSystem, Num} <: AbstractCoconVector
#  cs::C
#  Aᵢ::Vector{Num} # for a component Aᵢ uⁱ, this represents Aᵢ
#end
Base.getindex(a::CovariantVector, i) = a.Aⁱ[i]
Base.setindex!(a::CovariantVector, v, i) = (a.Aⁱ[i] = v)
Base.getindex(a::ContravariantVector, i) = a.Aᵢ[i]
Base.setindex!(a::ContravariantVector, v, i) = (a.Aᵢ[i] = v)

function Base.convert(::Type{ContravariantVector}, a::CovariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, gⁱʲ(cs) * a.Aᵢ)
end

function Base.convert(::Type{CovariantVector}, a::ContravariantVector)
  cs = coordinatesystem(a)
  return CovariantVector(cs, gᵢⱼ(cs) * a.Aⁱ)
end

Base.convert(::Type{ContravariantVector}) = CovariantVector
Base.convert(::Type{CovariantVector}) = ContravariantVector

function Base.:*(a::Num, b::ContravariantVector)
  return ContravariantVector(coordinatesystem(b), a .* b.Aⁱ)
end
Base.:*(a::AbstractCoconVector, b::Num) = b * a
for (CV, Ai) in ((CovariantVector, :Aᵢ), (ContravariantVector, :Aⁱ))
  for op in (:+, :-)
    @eval function Base.$(op)(a::T, b::T) where {T<:$(CV)}
      cs = coordinatesystem(a)
      @assert cs == coordinatesystem(b)
      return CovariantVector(cs, a.$(Ai) + b.$(Ai))
    end
  end
  for op in (:*, :/)
    @eval function Base.$(op)(a::$(CV), b::Num)
      return CovariantVector(coordinatesystem(b), op.(a.$(Ai), b))
    end
  end
end
Base.:*(a::Num, b::AbstractCoconVector) = b * a

coordinatesystem(a) = a.cs
gᵢⱼ(c::AbstractCoordinateSystem) = c.∂u⃗_∂x⃗ * c.∂u⃗_∂x⃗'
gⁱʲ(c::AbstractCoordinateSystem) = inv(gᵢⱼ(c))
𝐽(c::AbstractCoordinateSystem) = sqrt(det(gᵢⱼ(c)))

function (d::∇)(q)
  return ContravariantVector(d.c, [d.c.∂₁(q), d.c.∂₂(q), d.c.∂₃(q)])
end

function (d::∇o)(a::ContravariantVector)
  cs = coordinatesystem(a)
  J = 𝐽(cs)
  f = (cs.∂₁(J * a.Aⁱ[1]) + cs.∂₂(J * a.Aⁱ[2]) + cs.∂₃(J * a.Aⁱ[3])) / J
  return Symbolics.expand_derivatives(f)
end
∇o(a::CovariantVector) = ∇o(convert(a))

function (d::∇x)(a::CovariantVector)
  cs = coordinatesystem(a)
  J = 𝐽(cs)
  return ContravariantVector(cs, [cs.∂₂(a.Aᵢ[3]) - cs.∂₃(a.Aᵢ[2]) / J,
                                  cs.∂₃(a.Aᵢ[1]) - cs.∂₁(a.Aᵢ[3]) / J,
                                  cs.∂₁(a.Aᵢ[2]) - cs.∂₂(a.Aᵢ[1]) / J])
end
∇x(a::ContravariantVector) = ∇x(convert(CovariantVector, a))

import LinearAlgebra: cross, dot
dot(a::ContravariantVector, b::CovariantVector) = dot(a.Aⁱ, b.Aᵢ)
dot(a::CovariantVector, b::ContravariantVector) = dot(b, a)
dot(a::T, b::T) where {T<:AbstractCoconVector} = dot(a, convert(convert(T), b))
function cross(a::CovariantVector, b::CovariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, cross(a, b) ./ 𝐽(cs))
end
function cross(a::ContravariantVector, b::ContravariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, cross(a.Aⁱ, b.Aⁱ) .* 𝐽(cs))
end
cross(a::CovariantVector, b::ContravariantVector) = cross(a, convert(CovariantVector, b))
cross(a::ContravariantVector, b::CovariantVector) = cross(convert(CovariantVector, a), b)
export cross, dot
end # module

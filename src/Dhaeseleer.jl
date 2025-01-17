module Dhaeseleer

using Symbolics, LinearAlgebra

sqrtrule1 = @rule sqrt((~x)^2) => ~x
sqrtrule2 = @rule sqrt(1/(~x)^2) => 1/(~x)

function process(a)
  a = Symbolics.expand_derivatives(a, true)
  a = Symbolics.expand(a)
  a = Symbolics.simplify(a, Rewriters.Prewalk(Rewriters.PassThrough(sqrtrule1)))
  a = Symbolics.simplify(a, Rewriters.Prewalk(Rewriters.PassThrough(sqrtrule2)))
  return a
end

abstract type AbstractCoordinateSystem end

struct ArbitraryCoordinateSystem <: AbstractCoordinateSystem
  xi::Vector{Num} # Cartesian coordinates
  ui::Vector{Num} # Non-Cartesian coordinates
  ∂₁::Differential # Differential with respect to 1st non-Cartesian coord
  ∂₂::Differential # Differential with respect to 2nd non-Cartesian coord
  ∂₃::Differential # Differential with respect to 3rd non-Cartesian coord
  ∂u⃗_∂x⃗::Matrix{Num} # 3x3 matrix of ∂ui/∂xj where i (j) are rows (columns)
end

"""
    gⁱʲ(c::AbstractCoordinateSystem)

Return the metric matrix of (uⁱ)ᵀuʲ
"""
gⁱʲ(c::AbstractCoordinateSystem) = c.∂u⃗_∂x⃗ * c.∂u⃗_∂x⃗'
"""
    gᵢⱼ(c::AbstractCoordinateSystem)

Return the metric matrix of (uᵢ)ᵀuⱼ
"""
gᵢⱼ(c::AbstractCoordinateSystem) = inv(gⁱʲ(c))
"""
    jac(c::AbstractCoordinateSystem)

Return the scalar Jacobian of the coordinate system
"""
jac(c::AbstractCoordinateSystem) = sqrt(det(gᵢⱼ(c)))

abstract type AbstractDifferentialOperator end
for f in (:gⁱʲ, :gᵢⱼ, :jac)
  @eval $f(d::AbstractDifferentialOperator) = $f(coordinatesystem(d))
end

"""
    Gradient operator ∇

# Fields
- c: AbstractCoordinateSystem
"""
struct ∇ <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end

"""
    Divergence operator ∇o

# Fields
- c: AbstractCoordinateSystem
"""
struct ∇o <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end
"""
    Curl operator ∇x

# Fields
- c: AbstractCoordinateSystem
"""
struct ∇x <: AbstractDifferentialOperator
  c::AbstractCoordinateSystem
end
coordinatesystem(a::AbstractDifferentialOperator) = a.c


"""
    CoordinateSystem(xi,ui)

Pass in a vector of symbols for the underlying Cartesian coordinates `xi`
and the symbolic coordinates for another coordinate system `ui`.

...
# Arguments
- `xi`:  Cartesian coordinates
- `ui`:  A Vector of symbols for another coordinate system
...

# Example
```julia
@variables x y z R ϕ Z
xi = [x, y, z]
ui = [R, ϕ, Z]
∇, ∇o, ∇x = CoordinateSystem(xi, ui)
```
"""
function CoordinateSystem(xi, ui)
  ∂₁ = Differential(ui[1])
  ∂₂ = Differential(ui[2])
  ∂₃ = Differential(ui[3])
  ∂u⃗_∂x⃗ = Matrix{Num}(undef, 3, 3)
  ss(x) = String(Symbol(x))
  ∂x = Differential(xi[1])
  ∂y = Differential(xi[2])
  ∂z = Differential(xi[3])
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

"""
    subsimp!(input::AbstractArray,dict)

Substitute a symbols defined as pairs in a dictionary into an abstract
array and then simplify the expressions within afterwards.

This is an in-place operation and modifies the contents of the array.

...
# Arguments
- `input::AbstractArray`: 
- `dict`: 
...

```julia
using Dhaeseleer, Symbolics
@variables x::Real, y::Real, z::Real
@variables R ϕ Z
∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
∇ϕ = ∇(ϕ)
d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
Dhaeseleer.subsimp!(∇ϕ, d)
```
"""
function subsimp!(input::AbstractArray, dict)
  for i in eachindex(input)
    input[i] = Symbolics.substitute(input[i], dict)
  end
  return simplify!(input)
end

"""
    subsimp!(input::Num,dict)

Return a symbol after substituting symbols in it according to the contents
of a dictionary and simplify the result.

This is creates a copy of the incoming symbolic expression

...
# Arguments
- `input::Num`: 
- `dict`: 
...

```julia
using Dhaeseleer, Symbolics
@variables x::Real, y::Real
@variables R
Dhaeseleer.subsimp(sqrt(x^2 + y^2), x^2=>R^2 - y^2)
```
"""
subsimp(input::Num, dict) = simplify(Symbolics.substitute(input, dict))

coordinatesystem(a) = a.cs

"""
    Base.convert(a::CovariantVector)

Switcha CovariantVector to a ContravariantVector

...
# Arguments
- `a::CovariantVector`: 
...
"""
function Base.convert(a::CovariantVector)
  cs = coordinatesystem(a)
  return ContravariantVector(cs, gⁱʲ(cs) * a.Aᵢ)
end

"""
    Base.convert(a::ContravariantVector)

Switcha ContravariantVector to a CovariantVector

...
# Arguments
- `a::ContravariantVector`: 
...
"""
function Base.convert(a::ContravariantVector)
  cs = coordinatesystem(a)
  return CovariantVector(cs, gᵢⱼ(cs) * a.Aⁱ)
end

"""
    unitise(a::ContravariantVector)

Return a vector with unit vectors in the coordinate system.
The coefficients have physical units.

...
# Arguments
- `a::ContravariantVector`: 
...

```julia
using Dhaeseleer, Symbolics, Test
@variables x::Real, y::Real, z::Real
@variables R ϕ Z
∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
∇ϕ = ∇(ϕ)
d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
Aᵢuⁱ = Dhaeseleer.convert(∇ϕ)
Aϕ̂ = Dhaeseleer.unitise(Aᵢuⁱ)
Dhaeseleer.subsimp!(Aϕ̂, d)
Dhaeseleer.subsimp!(Aϕ̂, x^2=>R^2 - y^2)
@test all(Aϕ̂ .- [0, 1/R, 0] .== zeros(3))
```
"""
function unitise(a::ContravariantVector)
  cs = coordinatesystem(a)
  #Aⁱ uᵢ = Ai ei, ei = uᵢ / |uᵢ|, Aⁱ uᵢ = Ai uᵢ / |uᵢ|, hence Aⁱ |uᵢ| = Ai
  gij = gᵢⱼ(cs)
  return [a.Aⁱ[1] * sqrt(gij[1, 1]),
          a.Aⁱ[2] * sqrt(gij[2, 2]),
          a.Aⁱ[3] * sqrt(gij[3, 3])]
end

"""
    unitise(a::CovariantVector)

Return a vector with unit vectors in the coordinate system.
The coefficients have physical units.

...
# Arguments
- `a::CovariantVector`: 
...

```julia
using Dhaeseleer, Symbolics, Test
@variables x::Real, y::Real, z::Real
@variables R ϕ Z
∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
∇ϕ = ∇(ϕ)
d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
n∇ϕ = Dhaeseleer.unitise(∇ϕ)
Dhaeseleer.subsimp!(n∇ϕ, d)
Dhaeseleer.subsimp!(n∇ϕ, x^2=>R^2 - y^2)
@test all(n∇ϕ .- [0, 1/R, 0] .== zeros(3))
```
"""
function unitise(a::CovariantVector)
  cs = coordinatesystem(a)
  #Aᵢ uⁱ = Ai ei, ei = uⁱ / |uⁱ|, Aᵢ uⁱ = Ai uⁱ / |uⁱ|, hence Aᵢ |uⁱ| = Ai
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

(d::∇)(q) = CovariantVector(d.c, [d.c.∂₁(q), d.c.∂₂(q), d.c.∂₃(q)])

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

struct Ao∇ <: AbstractDifferentialOperator
  A::ContravariantVector
  d::∇
  function Ao∇(A::ContravariantVector, d::∇)
    @assert coordinatesystem(A) == coordinatesystem(d)
    return new(A, d)
  end
end

import LinearAlgebra: dot
dot(A::ContravariantVector, d::∇) = Ao∇(A, d)
dot(A::CovariantVector, d::∇) = Ao∇(convert(A), d)

"""
    (d::Ao∇)(a::Num)

Apply the directional derivative operator to a symbolic expression.

...
# Arguments
- `d::Ao∇(a::Num`: 
...

# Example
```julia
```
"""
function (d::Ao∇)(a::Num)
  return d.A.Aⁱ[1] * d.d.∂₁(a) + d.A.Aⁱ[2] * d.d.∂₂(a) + d.A.Aⁱ[3] * d.d.∂₃(a)
end

dot(a::ContravariantVector, b::CovariantVector) = dot(a.Aⁱ, b.Aᵢ)
dot(a::CovariantVector, b::ContravariantVector) = dot(b, a)
dot(a::T, b::T) where {T<:AbstractCoconVector} = dot(a, convert(b))
import LinearAlgebra: cross
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


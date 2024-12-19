using Dhaeseleer, Symbolics
using Test

@testset "Dhaeseleer.jl" begin
  @variables x::Real, y::Real, z::Real
  @variables R ϕ Z
  ∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
  @variables B₀ R₀ κ₀ q₀ τ₀ a₀ σ₀

  X = (R - R₀) / a₀
  Y = Z / R₀
  E₀ = a₀ / R₀
  Ψ = (X-E₀/2*(1-X^2))^2 + (1-E₀/4)*(1+τ₀*E₀*X*(2+E₀*X))*Y^2/σ₀^2

#  Ψ = B₀ / (2 * R₀^2 * κ₀ * q₀) * (R * Z + κ₀^2 / 4 * (R^2 - R₀^2))

  ∇ϕ = ∇(ϕ)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(∇ϕ, d)
  Dhaeseleer.subsimp!(∇ϕ, x^2=>R^2 - y^2)
  n∇ϕ = Dhaeseleer.unitise(∇ϕ)
  Dhaeseleer.subsimp!(n∇ϕ, d)
  Dhaeseleer.subsimp!(n∇ϕ, x^2=>R^2 - y^2)
  @test all(n∇ϕ .- [0, 1/R, 0] .== zeros(3))

  gij = Dhaeseleer.gⁱʲ(∇)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(gij, d)
  Dhaeseleer.subsimp!(gij, x^2=>R^2 - y^2)
  @test all(gij .- [1 0 0; 0 1/R^2 0; 0 0 1] .== zeros(3, 3))

  gij = Dhaeseleer.gᵢⱼ(∇)
  Dhaeseleer.subsimp!(gij, d)
  Dhaeseleer.subsimp!(gij, x^2=>R^2 - y^2)
  @test all(gij .- [1 0 0; 0 R^2 0; 0 0 1] .== zeros(3, 3))

  ∇Ψ = ∇(Ψ)
  ∇x∇Ψ = ∇x(∇Ψ)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(∇x∇Ψ, d)
  d = Dict(x^2=>R^2 - y^2, y/x=>tan(ϕ), z=>Z)
  Dhaeseleer.subsimp!(∇x∇Ψ, d)
  @test all(∇x∇Ψ .== [0, 0, 0])

  ∇ϕ = ∇(ϕ)
  ∇x∇ϕ = ∇x(∇ϕ)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(∇x∇ϕ, d)
  Dhaeseleer.subsimp!(∇x∇ϕ, x^2=>R^2 - y^2)
  @test all(∇x∇ϕ .== [0, 0, 0])

  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), R^2=>x^2 + y^2, Z=>z)
  ∇RZ = ∇(R * Z)
  t = ∇x(∇RZ)
  Dhaeseleer.subsimp!(t, d)
  Dhaeseleer.subsimp!(t, x^2=>R^2 - y^2)
  @test all(t .== [0, 0, 0])

  fΨ = B₀ * R₀
  B = fΨ * ∇(ϕ) + cross(∇(Ψ), ∇(ϕ))
  J = ∇x(B)
  ∇oJ = ∇o(J)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  ∇oJ = Dhaeseleer.subsimp(∇oJ, d)
  @test ∇oJ == 0

  T = cross(∇(Ψ), B) / abs(B)
  To∇Ψ = dot(T, ∇Ψ)
  To∇Ψ = Symbolics.substitute(To∇Ψ, [R=>sqrt(x^2 + y^2), ϕ=>atan(y, x), Z=>z])
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  To∇Ψ = Dhaeseleer.subsimp(To∇Ψ, d)
  To∇Ψ = Dhaeseleer.subsimp(To∇Ψ, x^2=>R^2 - y^2)
  @test To∇Ψ == 0

  s = dot(T, ∇x(T)) # shear
  c = cross(∇x(B / abs(B)), B / abs(B)) # curvature
  Bo∇ = dot(B, ∇)
end

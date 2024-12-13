using Dhaeseleer, Symbolics
using Test

@testset "Dhaeseleer.jl" begin
  @variables x, y, z
  R = sqrt(x^2 + y^2)
  ϕ = atan(y, x)
  Z = z
  ∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
  @variables B₀ R₀ κ₀ q₀

  Ψ = B₀ / (2 * R₀^2 * κ₀ * q₀) * (R * Z + κ₀^2 / 4 * (R^2 - R₀^2))
  ∇Ψ = ∇(Ψ)
  ∇ϕ = ∇(ϕ)
  fΨ = B₀ * R₀
  B = fΨ * ∇ϕ + cross(∇Ψ, ∇ϕ)
  J = ∇x(B)
  @test ∇o(J) == 0

end

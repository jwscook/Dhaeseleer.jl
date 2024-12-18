using Dhaeseleer, Symbolics
using Test

@testset "Dhaeseleer.jl" begin
  @variables x::Real, y::Real, z::Real
  @variables R(x, y,z)::Real ϕ(x, y, z)::Real Z(x, y, z)::Real
  #@variables R ϕ Z
  ∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
  @variables B₀ R₀ κ₀ q₀

  Ψ = B₀ / (2 * R₀^2 * κ₀ * q₀) * (R * Z + κ₀^2 / 4 * (R^2 - R₀^2))
  ∇Ψ = ∇(Ψ)
  ∇ϕ = ∇(ϕ)

  gij = Dhaeseleer.gⁱʲ(∇.c)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(gij, d)
  Dhaeseleer.subsimp!(gij, x^2=>R^2 - y^2)
  @test all(gij .- [1 0 0; 0 1/R^2 0; 0 0 1] .== zeros(3, 3))

  gij = Dhaeseleer.gᵢⱼ(∇.c)
  Dhaeseleer.subsimp!(gij, d)
  Dhaeseleer.subsimp!(gij, x^2=>R^2 - y^2)
  @test all(gij .- [1 0 0; 0 R^2 0; 0 0 1] .== zeros(3, 3))

  ∇x∇Ψ = ∇x(∇Ψ)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  Dhaeseleer.subsimp!(∇x∇Ψ, d)
  d = Dict(x^2=>R^2 - y^2, y/x=>tan(ϕ), z=>Z)
  Dhaeseleer.subsimp!(∇x∇Ψ, d)
  @test all(∇x∇Ψ .== [0, 0, 0])


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

end

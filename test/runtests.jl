using Dhaeseleer, Symbolics
using Test

@testset "Dhaeseleer.jl" begin
  @variables x::Real, y::Real, z::Real
  #R = sqrt(x^2 + y^2)
  #ϕ = atan(y, x)
  #Z = z
  #@variables R(x, y,z)::Real ϕ(x, y, z)::Real Z(x, y, z)::Real
  @variables R ϕ Z
  ∇, ∇o, ∇x = Dhaeseleer.CoordinateSystem([x, y, z], [R, ϕ, Z])
  @variables B₀ R₀ κ₀ q₀

  Ψ = B₀ / (2 * R₀^2 * κ₀ * q₀) * (R * Z + κ₀^2 / 4 * (R^2 - R₀^2))
  ∇Ψ = ∇(Ψ)
  ∇ϕ = ∇(ϕ)
  ∇x∇Ψ = ∇x(∇Ψ)
  Dhaeseleer.simplify(∇x∇Ψ)
  #@show ∇x∇Ψ.Aⁱ[1]
  #∇x∇Ψ.Aⁱ[1] = Symbolics.substitute(∇x∇Ψ.Aⁱ[1],
  #∇x∇Ψ.Aⁱ[1] = Symbolics.simplify(∇x∇Ψ.Aⁱ[1])
  #∇x∇Ψ.Aⁱ[1] = Symbolics.expand_derivatives(∇x∇Ψ.Aⁱ[1], true)
  #d = Dict(x^2=>R^2 - y^2, y/x=>tan(ϕ))
  #Dhaeseleer.subsimp!(∇x∇Ψ, d)
  #d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y, z), Z=>z)
  d = Dict(R=>sqrt(x^2 + y^2), ϕ=>atan(y / x), Z=>z, R^2=>x^2 + y^2)
  #d = Dict(R=>x, ϕ=>y, Z=>z)
  Dhaeseleer.subsimp!(∇x∇Ψ, d)
  #d = Dict((x^2 + y^2)=>R^2, y/x=>tan(ϕ))
  #Dhaeseleer.subsimp!(∇x∇Ψ, d)
  for _ in 1:1
    #d = Dict(sqrt(x^2 + y^2)=>R, y/x=>tan(ϕ))
    #Dhaeseleer.subsimp!(∇x∇Ψ, d)
    d = Dict(x^2=>R^2 - y^2, y/x=>tan(ϕ), z=>Z)
    #d = Dict(x=>R, y=>ϕ, z=>Z)
    Dhaeseleer.subsimp!(∇x∇Ψ, d)
  end
  @show ∇x∇Ψ.Aⁱ
  #@show ∇x∇Ψ.Aⁱ[3]
  #fΨ = B₀ * R₀
  #B = fΨ * ∇ϕ + cross(∇Ψ, ∇ϕ)
  #J = ∇x(B)
  #∇oJ = ∇o(J)
  #∇oJ = Symbolics.substitute(∇oJ, [R=>sqrt(x^2 + y^2), ϕ=>atan(y, x), Z=>z])
  #@show expand_derivatives(J)
  #@show Dhaeseleer.expand_derivatives(∇oJ)
  #@test ∇oJ == 0
  #T = cross(∇Ψ, B) / abs(B)
  #To∇Ψ = dot(T, ∇Ψ)
  ##To∇Ψ = Symbolics.substitute(To∇Ψ, [R=>sqrt(x^2 + y^2), ϕ=>atan(y, x), Z=>z])
  #@show Dhaeseleer.expand_derivatives(To∇Ψ)
  #@test To∇Ψ == 0

end

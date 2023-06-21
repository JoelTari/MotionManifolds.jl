using ManifoldExtras
using Test

@testset "ManifoldExtras.jl" begin
    # Write your tests here.
    println("Declaring variables for tests")
    x=-2.3;y=6.5;th=-5*pi/6;th2=43*pi/6
    R=[cos(th) -sin(th); sin(th) cos(th)]
    # X = makeSE2(x,y,a)
    rand1=rand()
    @test R ≈ to_matrix(SO2(th))
    @test to_matrix(SO2(th+th2)) ≈ to_matrix(SO2(th)*SO2(th2))
    @test SO2(th+th2).th ≈ (SO2(th)*SO2(th2)).th
    @test SO2(th+th2).th != (SO2(th)*SO2(th2+2*pi)).th # th2 not the same +2pi
    @test SO2(th+th2).c ≈ (SO2(th)*SO2(th2+2*pi)).c
    @test ecpi(SO2(th)) ≈ ecpi(SO2(th-234*2*pi))
    @test inv(SO2(th2)).c ≈ cos(-th2)
    @test inv(SO2(th2)).s ≈ sin(-th2)
    @test Log(Exp(5.001*pi)) ≈ (5.001-6)*pi # non-bijectivity
    @test Log(Exp(rand1*8*pi)) ≈ ecpi(Exp(rand1*8*pi))
    # @test [cos(a) -sin(a); sin(a) cos(a)] ≈ getRotationComponent(X)        
end

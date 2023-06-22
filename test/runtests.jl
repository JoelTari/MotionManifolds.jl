using ManifoldExtras, StaticArrays
using Test

println("Declaring variables for tests")
x=-2.3;y=6.5;th=-5*pi/6;th2=43*pi/6
X = SE2(SA_F64[x,y],SO2(th))
Xb = SE2([x,y],SO2(th))
Xc = SE2([x,y],th)
# sk = log_lie(X)
R=[cos(th) -sin(th); sin(th) cos(th)]
tau = SA_F64[23.,-8,0.9]
tauOOB = SA_F64[23.,-8,0.9-4*pi]
taud = [tau...] # non static Vector{Float64}
tau0 = SA_F64[0,0,0]
rand1=rand()
sk = se2(taud)
sk = se2(tau)
Xsk = ManifoldExtras.exp_lie(sk)
# 
# to_matrix(Exp(tau))
# to_matrix(se2(tau))
# TODO: hat vee
# TODO: r/l jacobians

println(":: Starting the test set ::")

@testset "ManifoldExtras.jl" begin
    # Write your tests here.
    # X = makeSE2(x,y,a)
    @test R ≈ to_matrix(SO2(th))
    @test to_matrix(SO2(th+th2)) ≈ to_matrix(SO2(th)*SO2(th2))
    @test SO2(th+th2).th ≈ (SO2(th)*SO2(th2)).th
    @test SO2(th+th2).th != (SO2(th)*SO2(th2+2*pi)).th # th2 not the same +2pi
    @test SO2(th+th2).c ≈ (SO2(th)*SO2(th2+2*pi)).c
    @test ecpi(SO2(th)) ≈ ecpi(SO2(th-234*2*pi))
    # @test Exp(th) ≈ Exp(th+12pi) # should work for SO2 
    @test Log(Exp(th)) ≈ Log(Exp(th+12pi))
    @test inv(SO2(th2)).c ≈ cos(-th2)
    @test inv(SO2(th2)).s ≈ sin(-th2)
    @test Log(Exp(5.001*pi)) ≈ (5.001-6)*pi # non-bijectivity
    @test Log(Exp(rand1*8*pi)) ≈ ecpi(Exp(rand1*8*pi))
    @test Log(Exp(tau)) ≈ tau 
    @test Log(Exp(taud)) ≈ taud # note here: dynamic input vector
    @test Log(Exp(tauOOB)) != tauOOB
    # @test Exp(tau) ≈  
    @test Jr(se2(tau)) ≈ inv(Jrinv(se2(tau))) 
    @test Jr(se2(tau)) ≈ ExpAndJr(tau)[2]
    @test Jl(se2(tauOOB)) ≈ inv(Jlinv(se2(tauOOB))) 
    # @test_ X.rot = pi/6 # immutable
end

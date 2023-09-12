# Copyright 2023 AKKA Technologies and LAAS-CNRS (joel.tari@akkodis.com) 
# 
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by 
# the European Commission - subsequent versions of the EUPL (the "Licence"); 
# You may not use this work except in compliance with the Licence. 
# You may obtain a copy of the Licence at: 
# 
# https://joinup.ec.europa.eu/software/page/eupl 
# 
# Unless required by applicable law or agreed to in writing, software 
# distributed under the Licence is distributed on an "AS IS" basis, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and 
# limitations under the Licence. 

using ManifoldExtras, StaticArrays
using Test

println("Declaring variables for tests")
x = -2.3;
y = 6.5;
th = -5 * pi / 6;
th2 = 43 * pi / 6;
X = SE2(SA_F64[x, y], SO2(th))
Xb = SE2([x, y], SO2(th))
Xc = SE2([x, y], th)
# sk = log_lie(X)
R = [cos(th) -sin(th); sin(th) cos(th)]
tau = SA_F64[23.0, -8, 0.9]
tau0a = SA_F64[23.0, -8, 0]
tauOOB = SA_F64[23.0, -8, 0.9-4*pi]
taud = [tau...] # non static Vector{Float64}
tau0 = SA_F64[0, 0, 0]
rand1 = rand()
sk = se2(taud)
sk = se2(tau)
Xsk = ManifoldExtras.exp_lie(sk)[1]

println(":: Starting the test set ::")

@testset "ManifoldExtras.jl" begin
    @test R ≈ to_matrix(SO2(th))
    @test to_matrix(SO2(th + th2)) ≈ to_matrix(SO2(th) * SO2(th2))
    @test SO2(th + th2).th ≈ (SO2(th) * SO2(th2)).th
    @test SO2(th + th2).th != (SO2(th) * SO2(th2 + 2 * pi)).th # th2 not the same +2pi
    @test SO2(th + th2).c ≈ (SO2(th) * SO2(th2 + 2 * pi)).c
    @test ecpi(SO2(th)) ≈ ecpi(SO2(th - 234 * 2 * pi))
    @test Log(Exp(th)) ≈ Log(Exp(th + 12pi))
    @test inv(SO2(th2)).c ≈ cos(-th2)
    @test inv(SO2(th2)).s ≈ sin(-th2)
    @test Log(Exp(5.001 * pi)) ≈ (5.001 - 6) * pi # non-bijectivity
    @test Log(Exp(rand1 * 8 * pi)) ≈ ecpi(Exp(rand1 * 8 * pi))
    @test Log(Exp(tau)) ≈ tau
    @test Log(Exp(tau0)) ≈ tau0
    @test Log(Exp(tau0a)) ≈ tau0a
    @test Log(Exp(taud)) ≈ taud # note here: dynamic input vector
    @test Log(Exp(tauOOB)) != tauOOB
    @test Jr(se2(tau)) ≈ inv(Jrinv(se2(tau)))
    @test Jr(se2(tau)) ≈ ExpAndJr(tau)[2]
    @test Jl(se2(tauOOB)) ≈ inv(Jlinv(se2(tauOOB)))
    @test to_matrix(Exp(tau)) ≈ to_matrix(exp_lie(se2(tau))[1])
    @test to_matrix(hat(tau)) ≈ to_matrix(se2(tau))
    @test to_matrix(inv(X)) ≈ inv(to_matrix(X))
    @test to_matrix(inv(Xsk)) ≈ inv(to_matrix(Xsk))
    @test to_matrix(Xsk * inv(X)) ≈ to_matrix(Xsk) * inv(to_matrix(X))
    @test Jlinv(sk) ≈ inv(Jl(sk))
    @test Jlinv(se2(tau0)) ≈ inv(Jl(se2(tau0)))
    @test Jlinv(se2(tau0a)) ≈ inv(Jl(se2(tau0a)))
    @test Jrinv(sk) ≈ inv(Jr(sk))
    @test Jrinv(se2(tau0)) ≈ inv(Jr(se2(tau0)))
    @test Jrinv(se2(tau0a)) ≈ inv(Jr(se2(tau0a)))
    @test Jr(se2(tau0)) ≈ Jl(se2(-tau0))
    @test Jr(se2(tau0a)) ≈ Jl(se2(-tau0a))
    @test Jr(se2(tau)) ≈ Jl(se2(-tau))
    @test Jr(se2(tauOOB)) ≈ Jl(se2(-tauOOB))
    @test ExpAndJr(vee(sk))[2] ≈ Jr(sk)
    @test ExpAndJr(tau0a)[2] ≈ Jr(hat(tau0a))
    # @test_ X.rot = pi/6 # immutable
end

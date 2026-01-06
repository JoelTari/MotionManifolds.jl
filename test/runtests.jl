# BSD 2-Clause License
#
# Copyright (c) 2024, LAAS-CNRS & AKKODIS
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using MotionManifolds, StaticArrays
using Test

println("Declaring variables for tests")
x = -2.3;
y = 6.5;
th = -5 * pi / 6;
th2 = 43 * pi / 6;
X = SE2(SA_F64[x, y], SO2(th))
# Xc = SE2([x, y], th)
# sk = log_lie(X)
R = SMatrix{2,2,Float64,4}([cos(th) -sin(th); sin(th) cos(th)])
tau = SA_F64[23.0, -8, 0.9]
tau0a = SA_F64[23.0, -8, 0]
tauOOB = SA_F64[23.0, -8, 0.9-4*pi]
tau0 = SA_F64[0, 0, 0]
rand1 = rand()
sk = se2(tau)
Xsk = MotionManifolds.exp_lie(sk)[1]

println(":: Starting the test set ::")

@testset "MotionManifolds.jl" begin
    @test R ≈ to_matrix(SO2(th))
    @test to_matrix(SO2(R)) ≈ to_matrix(SO2(th))
    @test to_matrix(SO2(th + th2)) ≈ to_matrix(SO2(th) * SO2(th2))
    @test SO2(th + th2).th ≈ (SO2(th) * SO2(th2)).th
    @test SO2(th + th2).th != (SO2(th) * SO2(th2 + 2 * pi)).th # th2 not the same +2pi
    @test SO2(th + th2).c ≈ (SO2(th) * SO2(th2 + 2 * pi)).c
    @test ecpi(SO2(th)) ≈ ecpi(SO2(th - 234 * 2 * pi))
    @test Log(Exp(th,so2)) ≈ Log(Exp(th + 12pi,so2))
    @test inv(SO2(th2)).c ≈ cos(-th2)
    @test inv(SO2(th2)).s ≈ sin(-th2)
    @test Log(Exp(5.001 * pi,so2)) ≈ (5.001 - 6) * pi # non-bijectivity
    @test Log(Exp(rand1 * 8 * pi,so2)) ≈ ecpi(Exp(rand1 * 8 * pi,so2))
    @test Log(Exp(tau,se2)) ≈ tau
    @test Log(Exp(tau0,se2)) ≈ tau0
    @test Log(Exp(tau0a,se2)) ≈ tau0a
    @test Log(Exp(tauOOB,se2)) != tauOOB
    @test Jr(se2(tau)) ≈ inv(Jrinv(se2(tau)))
    # @test Jr(se2(tau)) ≈ ExpAndJr(tau)[2]
    @test Jl(se2(tauOOB)) ≈ inv(Jlinv(se2(tauOOB)))
    @test to_matrix(Exp(se2(tau))) ≈ to_matrix(exp_lie(se2(tau))[1])
    @test to_matrix(hat(tau,se2)) ≈ to_matrix(se2(tau))
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
    # @test ExpAndJr(vee(sk))[2] ≈ Jr(sk)
    # @test ExpAndJr(tau0a)[2] ≈ Jr(hat(tau0a))
    # @test_ X.rot = pi/6 # immutable

    # test Exp via MotionManifolds.numExp
    @test isapprox(MotionManifolds.numExp(Log(X),se2,18), to_matrix(X))

    # other test: the Adjoint
    res1=to_matrix(X)*to_matrix(se2(tau))*to_matrix(inv(X))
    res2=hat(Adjm(X)*vee(se2(tau)),se2) |> to_matrix
    @test isapprox(res1,res2)
    rot=SO2(th)
    res1=to_matrix(rot)*to_matrix(so2(th2))*to_matrix(inv(rot))
    res2=hat(Adjm(rot)*vee(so2(th2)),so2) |> to_matrix
    @test isapprox(res1,res2)

    stau=se2(tau)
    sth=so2(th)
    vr=SVector(rand(6)...)
    svr=se3(vr)
    xr=SVector(rand(3)...)
    sxr=so3(xr)
    @test to_matrix(Exp(stau)) ≈ to_matrix(Exp(tau, se2))
    @test to_matrix(Exp(sth)) ≈ to_matrix(Exp(th, so2))
    @test to_matrix(Exp(sxr)) ≈ to_matrix(Exp(xr, so3))
    @test to_matrix(Exp(svr)) ≈ to_matrix(Exp(vr, se3))

    @test Jlinv(svr) ≈ inv(Jl(svr))
    @test Jrinv(svr) ≈ inv(Jr(svr))
    @test Jlinv(sxr) ≈ inv(Jl(sxr))
    @test Jrinv(sxr) ≈ inv(Jr(sxr))
    # adjoint relations for SE3 and SO3
    Xvr=rand(SE3)
    res1=to_matrix(Xvr)*to_matrix(svr)*to_matrix(inv(Xvr))
    res2=hat(Adjm(Xvr)*vee(svr),se3) |> to_matrix
    Xr=rand(SO3)
    res1=to_matrix(Xr)*to_matrix(sxr)*to_matrix(inv(Xr))
    res2=hat(Adjm(Xr)*vee(sxr),so3) |> to_matrix
    @test isapprox(res1,res2)
end

@testset "MotionManifolds operators" begin
    # test + * dot operators
    th1,th2=rand(2)
    @test SO2(th1)*SO2(th2) == SO2(th1)+SO2(th2)
    X1,X2=[SE2(5*randn(2)...,th1),SE2(2*randn(2)...,th2)]
    @test X1+X2 == X1*X2
    Xr1, Xr2 = SO3(SA_F64[3*randn(3)...]),SO3(SA_F64[randn(3)...])
    @test Xr1+Xr2 == Xr1*Xr2
    @test Xr1+Xr2 != Xr2+Xr1 # no commutativity
    Y1,Y2 = SE3(SA_F64[5*randn(3)...],Xr1), SE3(SA_F64[randn(3)...],Xr2)
    @test Y1+Y2 == Y1*Y2
    @test Y1+Y2 != Y2+Y1
    # minus SO2
    matth1minusth2 = to_matrix(SO2(th1)-SO2(th2))
    matth1minusth2_ = to_matrix(th1|>SO2)*to_matrix(inv(th2|>SO2))
    @test isapprox(matth1minusth2,matth1minusth2_)
    # minus SE2
    matX1minusX2 = to_matrix(X1-X2)
    matX1minusX2_ = to_matrix(X1)*to_matrix(inv(X2))
    @test isapprox(matX1minusX2,matX1minusX2_)
    @test X1-X2 != X2-X1
    # minus SO3
    matXr1minusXr2 = to_matrix(Xr1-Xr2)
    matXr1minusXr2_ = to_matrix(Xr1)*to_matrix(inv(Xr2))
    @test isapprox(matXr1minusXr2,matXr1minusXr2_)
    @test Xr1-Xr2 != Xr2-Xr1
    # minus SE3
    matY1minusY2 = to_matrix(Y1-Y2)
    matY1minusY2_ = to_matrix(Y1)*to_matrix(inv(Y2))
    @test isapprox(matY1minusY2,matY1minusY2_)
    @test Y1-Y2 != Y2-Y1
end

@testset "Transformation to/from Quaternions" begin
  # Tests for quaternions:
  #
  xw=so3(SA_F64[0.7937450088478044, 0.5582832527134615, 0.24143046756545883], 1.0282915415168719)
  W=Exp(vee(xw),SO3)
  # text 1 : check that q == q_bis (approx)
  q=to_quat(xw)
  q_bis=to_quat(W)
  @test q ≈ q_bis
  # text 2 : check those 4 matrices are the same (approx)
  @test to_matrix(q) ≈ to_matrix(W) ≈ to_matrix(SO3(q)) ≈ to_matrix(Exp(so3(q) |> vee, SO3))
  # quat operations
  q1=rand(Quaternion)
  q2=rand(Quaternion)
  W1 = SO3(q1)
  W2 = SO3(q2)
  @test isapprox(q1+q2 |> SO3 |> to_matrix, W1+W2 |> to_matrix) 
  @test isapprox(q1-q2 |> SO3 |> to_matrix, W1-W2 |> to_matrix) 
  W1=rand(SO3)
  W2=rand(SO3)
  q1=to_quat(W1)
  q2=to_quat(W2)
  @test isapprox(q1+q2 |> SO3 |> to_matrix, W1+W2 |> to_matrix) 
  @test isapprox(q1-q2 |> SO3 |> to_matrix, W1-W2 |> to_matrix) 
  # quaternion action on a 3D vector
  W=rand(SO3)
  q=to_quat(W)
  x=rand(SVector{3,Float64})
  @test isapprox( W+x, q+x )
end

@testset "SO3 Matrix" begin
  # test against a former bug (encountered in real life)
  R= SA_F64[-0.999764194980216 0.0012702695881678077 -0.021678119169681133;
           0.0012700454281182416 -0.9931583301449304 -0.11676865267859228;
          -0.021678132303594052 -0.11676865024027544 0.9929225251251368]
  @test isapprox( to_matrix(SO3(R)), R)
  R=SA_F64[-1 0 0.0;0 -1 0;0 0 1]
  @test isapprox( to_matrix(SO3(R)), R)
  R=SA_F64[1.0 0 0;0 1 0;0 0 1]
  @test isapprox( to_matrix(SO3(R)), R)
  uw=SA_F64[0.611,0.532,0.962]
  R=SO3(uw).R
  Xr=SO3(R)
  @test isapprox( Xr.u*Xr.w, uw)
end

@testset "Manifold - (minus associated methods)" begin
  Z=SE2(rand(3)...)
  Z2=SE2(rand(3)...)
  @test -Z == inv(Z)
  @test isapprox(-Z+Z2 |> to_matrix , inv(Z)*Z2 |> to_matrix)
  @test isapprox(Z-Z2 |> to_matrix, Z*inv(Z2) |> to_matrix)
  # SE3
  Zo=SO3(SVector(rand(3)...))
  Z=SE3(SVector(rand(3)...), Zo)
  Zo2=SO3(SVector(rand(3)...))
  Z2=SE3(SVector(rand(3)...), Zo2)
  @test -Z == inv(Z)
  @test isapprox(-Z+Z2 |> to_matrix , inv(Z)*Z2 |> to_matrix)
  @test isapprox(Z-Z2 |> to_matrix , Z*inv(Z2) |> to_matrix)
end

@testset "Lie algebra - (minus associated methods)" begin
  Z=SE2(rand(3)...)
  zeta=se2(SVector(rand(3)...))
  @test Z+Exp(vee(-zeta), se2) == Z+Exp(-vee(zeta), se2) == Z-zeta
  #
  Zo=SO3(SVector(rand(3)...))
  Z3=SE3(SVector(rand(3)...), Zo)
  zeta=se3(SVector(rand(6)...))
  @test Z3+Exp(vee(-zeta), se3) == Z3+Exp(-vee(zeta),se3) == Z3-zeta
end

@testset "Adjm for trivial groups and zero constructors" begin
  @test Adjm(2.0)==1
  @test Adjm(SA_F64[1.,1.,1.]) - [1 0 0;0 1 0;0 0 1] == zeros(3,3)
  # @test AdjmZero(SVector{2, Float64}) == Adjm(SVector{2,Float64}(rand(2)))
  @test AdjmZero(Float64) == 1
  @test AdjmZero(SVector{2, Float64}) - [1 0;0 1] == zeros(2,2)
  @test AdjmZero(SE2) - [1 0 0;0 1 0;0 0 1] == zeros(3,3)
  # Adjm(SVector{2,Float64}(rand(2)))
end

@testset "action" begin
  M=rand(SE2); p=rand(SVector{2})
  @test M+p ≈  M*p
  M=rand(SO2); p=rand(SVector{2})
  @test M+p ≈  M*p
  M=rand(SO3); p=rand(SVector{3})
  @test M+p ≈ M*p
  M=rand(SE3); p=rand(SVector{3})
  @test M+p ≈ M*p
end

@testset "isLieAlgebra" begin
 @test isLieAlgebra(se2,SE2)
 @test !isLieAlgebra(se2,Float64)
 @test !isLieAlgebra(se2,SE3)
 @test isLieAlgebra(se3,SE3)
 @test isLieAlgebra(so3,SO3)
 @test !isLieAlgebra(so2,SO3)
 @test isLieAlgebra(so2,SO2)
 @test isLieAlgebra(Float64,Float64)
 @test isLieAlgebra(SVector{2,Float64},SVector{2,Float64})
 @test !isLieAlgebra(SVector{2,Float64},SVector{3,Float64})
end

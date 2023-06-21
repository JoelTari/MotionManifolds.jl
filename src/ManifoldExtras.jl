module ManifoldExtras

# using Manifolds
# using RecursiveArrayTools
using StaticArrays
# add another dependency to confirm the workflow
# RecursiveArrayTools
# StaticArrays

export SO2, so2, SE2, se2, to_matrix, hat, vee, Adjm, ecpi, Log, Exp

# """
#     greet()
#
# Return a string that greets people.
#
# # Example
# ```julia
# julia> greet()
# "Hello ManifoldExtras."
# ```
#
# # To find out more
# See also [documentation](https://docs.julialang.org/en/v1/manual/documentation/).
# """
# function greet()
#   "Hello ManifoldExtras."
# end
#
# "The Special Euclidian 2 Manifold: \$\\mathrm{SE}(n)\$"
# const SE2LieGroup=SpecialEuclidean(2)
# const SE3LieGroup=SpecialEuclidean(2)
# const SO2LieGroup=SpecialOrthogonal(2)
# const SO3LieGroup=SpecialOrthogonal(2)
# # const fkso2=SMatrix{2,2,Float64,4}
# # const fkse2=SMatrix{3,3,Float64,9}
#
# """
#     NV
#
# NV is the number of nodes ``n``, valued 6.
#
# # Example
# ```jldoctest; setup = :(using ManifoldExtras)
# julia> NV
# 6
# ```
# """
# const NV = 6
#
# const is_point = Manifolds.is_point # TODO: export Manifolds instead
#
# "Type alias"
# # const SE2=ArrayPartition{Float64, Tuple{SVector{2, Float64}, SMatrix{2,2,Float64, 4}}}
# # const SE3=ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}}
# # const SO2=SMatrix{2,2,Float64,4}
# # const SO3=SMatrix{3,3,Float64,9}

struct SO2
  th::Float64 # not necessary inside [-pi,pi]
  c::Float64 # cos
  s::Float64 # sin
end
SO2(th::Float64) = SO2(th,cos(th),sin(th))
SO2(R::SMatrix{2,2,Float64}) = SO2(atan(R[2,1]/R[1,1]),R[1,1],R[2,1])

struct SE2
  t::SVector{2,Float64}
  rot::SO2
end
SE2(x::Float64,y::Float64,th::Float64) = SE2([x,y],SO2(th))

struct so2
  w::Float64
end

struct se2
  # a.k.a. screw
  vx::Float64
  vy::Float64
  w::Float64
end

ecpi(r::SO2) = atan(r.s,r.c) 
# force representation between [-pi,pi]; r.th may be outside (use case: external readability)

to_matrix(r::SO2) = SA_F64[r.c -r.s;r.s r.c]
to_matrix(X::SE2) = SA_F64[to_matrix(X.rot) X.t;0 0 1]
to_matrix(w::so2) = SA_F64[0 -w.w;w.w 0]
to_matrix(sk::se2) = SA_F64[to_matrix(sk.w) [sk.vx;sk.vy];0 0 0]

vee(w::so2) = w.w
vee(sk::se2) = SVector{3,Float64}[sk.vx sk.vy sk.w]
hat(w::Float64) = so2(w)
hat(tau::SVector{3,Float64}) = se2(tau...)

import Base: inv
Base.:inv(rot::SO2) = SO2(-rot.th,rot.c,-rot.s)
Base.:inv(X::SE2) = SE2(-inv(X.rot)*X.t,inv(X.rot))

import Base: *
function Base.:*(rot1::SO2,rot2::SO2)  
  SO2(rot1.th+rot2.th)
end

function Base.:*(X1::SE2,X2::SE2)
  SE2(X1.rot*X2.rot,X1.t + X1.rot*X2.t)
end


Adjm(rot::SO2) = 1

function Adjm(X::SE2)
SA_F64[X.rot.c -X.rot.s   -X.t[2]
       X.rot.s  X.rot.c    X.t[1]
       0        0          1      ]
end

struct Jexp end

exp_lie(w::so2)=SO2(w.w)
Exp(w::Float64) = exp_lie(hat(w))


function exp_lie(sk::se2)
  # technique from manifcpp
  # credit: sola/deray
  w_sq = sk.w*sk.w    
  # K1 is  sin(w)÷w 
  # K2 is  (1-cos(w))÷w     
  if w_sq < eps(Float64)
    K1 = 1 - w_sq/6
    K2 = .5*sk.w - 1.0/24*sk.w*w_sq
  else
    K1 = sin(sk.w)/sk.w;
    K2 = (1-cos(sk.w))/sk.w;
  end
  t = [K1 -K2;K2 K1]*[sk.vx;sk.vy];
  SE2(t,Exp(sk.w)), K1, K2, w_sq
end

Exp(tau::SVector{3,Float64}) = exp_lie(hat(tau))[1] 

function Exp(tau::SVector{3,Float64}, Nothing)
  expmap_value, K1 , K2, w_sq = exp_lie(hat(tau))
  vx = tau[1]; vy=tau[2]; w=tau[3]
  if w_sq < eps(Float64)
    J13 = -vy/2.0+w*vx/6.0
    J23 =  vx/2.0+w*vy/6.0
  else
    sw = sin(w)
    cw = cos(w)
    J13 = (-vy+w*vx+vy*cw-vx*sw)/w_sq
    J23 = ( vx+w*vy-vx*cw-vy*sw)/w_sq
  end
  Jr_exp = SA_F64[K1 K2 J13;-K2 K1 J23;0 0 1]
  expmap_value, Jr_exp
end

log_lie(rot::SO2) = so2(atan(rot.s,rot.c)) # non-bijectivity
Log(rot::SO2)=vee(log_lie(rot))


# TODO: Adjm SE2 SO2
# TODO: Exp Log (with optional jacobian)
#
# TODO: to_quat for SO3, to_S1 (unit circle) for SO2

# #multiply
#
# "Construct a SO2 object, compatible with Manifolds.jl"
# function makeSO2(a::Number)
#   aa=Float64(a)
#   SA_F64[cos(aa) -sin(aa);sin(aa) cos(aa)]
# end
#
#
# "get angle from an SO2 object, compatible with Manifolds.jl"
# function angleSO2(R::Union{Matrix{Float64},SMatrix{2,2,Float64}})
#   atan(R[2,1],R[1,1])
# end
#
# "extract translation vector from SE2"
# function getTranslationComponent(X::SE2)
#   SA_F64[X[1:2]...]
# end
#
# "extract rotation matrix from SE2"
# function getRotationComponent(X::SE2)
#   # reshape(SA_F64[X[3:end]...],2,2)
#   SMatrix{2,2,Float64,4}([X[3:end]...])
# end
#
# # "extract translation vector from SE3"
# # function getTranslationComponent(X::ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}})
# # end
# # "extract rotation matrix from SE3"
# # function getRotationComponent(X::ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}})
# #   reshape(SA_F64[X[3:end]...],2,2)
# # end
#
#
# "get angle from an SE2 object, compatible with Manifolds.jl"
# function angleSE2(X::SE2)
#   angleSO2(getRotationComponent(X))
# end
#
#
# "Construct an SE2 object from x,y,a."
# function makeSE2(x::Number,y::Number,a::Number)
#   ArrayPartition(SA_F64[Float64(x),Float64(y)],makeSO2(Float64(a)))
# end
#
# "Construct an SE2 object from (t, angle), compatible with Manifolds.jl"
# function makeSE2(t,a::Number)
#   ArrayPartition(SA_F64[t...],makeSO2(Float64(a)))
# end
#
# "Construct an SE2 object from (t, R), compatible with Manifolds.jl"
# function makeSE2(t,R::SO2)
#   ArrayPartition(SA_F64[t...],R)
# end
#
# "hat operation"
# function hat(a::Float64)
#   SA_F64[0 -a;a 0]
# end
# "vee operation"
# function vee(Sk::SMatrix{2,2,Float64,4})
#   Sk[2,1]
# end
#
# # "hat operation"
# # function hat(a::Float64)
# #   SA_F64[0 -a;a 0]
# # end
# # "vee operation"
# # function vee(Sk::SMatrix{2,2,Float64,4})
# #   Sk[2,1]
# # end
#
# # function 
# #
# # end
# #
# # hat
# # vee
# # Exp
# # Log
# # Adjm adjoint matrix
# #
# # rpySE3 quatSE3

end

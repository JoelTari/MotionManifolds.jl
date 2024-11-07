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

module MotionManifolds

using StaticArrays

export SO2,
    so2,
    SE2,
    se2,
    SO3,
    so3,
    SE3,
    se3,
    to_matrix,
    hat,
    wedge,
    vee,
    Adjm,
    ecpi,
    Log,
    Exp,
    log_lie,
    exp_lie,
    Jr,
    Jrinv,
    Jl,
    Jlinv,
    ExpAndJr,
    skew,
    is_skew

# export SO2FromMat

# julia> skew(SA_F64[1,2,3])
# 3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
#   0.0  -3.0   2.0
#   3.0   0.0  -1.0
#  -2.0   1.0   0.0

"""
    skew(a::SVector{3,Float64})

Compute the 3x3 skew matrix of a 3d vector

# Example
```jldoctest
julia> skew()
2×2 SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 0.0  -1.0
 1.0   0.0
```
"""
function skew(a::SVector{3,Float64})
  return SMatrix{3,3,Float64,9}(
      [  0    -a[3]  a[2]
        a[3]     0  -a[1]
       -a[2]   a[1]    0 ] 
      )
end

"""
    skew(a::Number=1)

Compute the 2x2 skew of a number (1 by default)
"""
function skew(a::Number=1)
  return SMatrix{2,2,Float64,4}(
    a*[0 -1;1 0]
  )
end

"""
    is_skew(M)

Check if a matrix is skewed.
"""
function is_skew(M)
  M'==-M
end


#--------------------------------------------------------------------#
#                                SO3                                 #
#--------------------------------------------------------------------#
# Possible constructors: SO3(); SO3(u,w,R); SO3(R); SO3(u,w); SO3(uw)
"""
    SO3

Special Orthogonal group 3

# Fields:
* `u::SVector{3,Float64}`: normed vector (direction)
* `w::Float64`: amplitude
* `R::SMatrix{3,3,Float64,9}`: Orthogonal Matrix.  ``det(R)=1``
"""
struct SO3
  u::SVector{3,Float64}
  w::Float64
  R::SMatrix{3,3,Float64,9}

  @doc """
      SO3()
  """
  function SO3()
    new(SVector{3,Float64}([1,0,0]),0,SMatrix{3,3,Float64,9}([1 0 0;0 1 0;0 0 1]))
  end
  #
  @doc """
      SO3(u,w,R)
  """
  function SO3(u::SVector{3,Float64},w::Float64,R::SMatrix{3,3,Float64,9})
    new(u,w,R)
  end
  #
  @doc """
      SO3(R)
  """
  function SO3(R::SMatrix{3,3,Float64,9})
    Rskew=R-R'
    @assert is_skew(Rskew)
    w=acos((R[1,1]+R[2,2]+R[3,3]-1)/2)
    if w > eps(Float64)
      uw=0.5*w/sin(w)*SA_F64[Rskew[3,2],Rskew[1,3],Rskew[2,1]]
      return new(uw/w,w,R)
    else
      uw=0.5*SA_F64[Rskew[3,2],Rskew[1,3],Rskew[2,1]]
      return new(uw,w,R)
    end
  end
  #
  @doc """
      SO3(u,w)
  """
  function SO3(u::SVector{3,Float64},w::Float64)
    R=[1 0 0;0 1 0;0 0 1]+sin(w)*skew(u)+(1-cos(w))*skew(u)^2
    new(u,w,R)
  end
  #
  @doc """
      SO3(uw)
  """
  function SO3(uw::SVector{3,Float64})
    w=sqrt(wu'wu)
    if w > eps(Float64)
      u=wu/sqrt(wu'wu)
    else
      u=wu
    end
    R=[1 0 0;0 1 0;0 0 1]+sin(w)*skew(u)+(1-cos(w))*skew(u)^2
    new(u,w,R)
  end
end


"""
    so3

``\\mathfrak{so}(3)``, the Lie algebra of the special orthogonal group 3 (SO3). Also known as rotational twist.

# Fields:
* `u::SVector{3,Float64}`: normed vector (direction)
* `w::Float64`: amplitude
"""
struct so3
  u::SVector{3,Float64}
  w::Float64

  @doc """
      so3()
  """
  function so3()
    new(SVector{3,Float64}([1,0,0]),0)
  end
  @doc """
      so3(uw)
  """
  function so3(uw::SVector{3,Float64})
    w=sqrt(uw'uw)
    w > eps(Float64) ? so3(uw/w,w) : so3(uw,w)
  end
end

"""
		vee(xr::so3)
"""
vee(xr::so3) = xr.u*xr.w
"""
		hat(uw::SVector{3,Float64}, ::Type{so3})::so3
"""
function hat(uw::SVector{3,Float64},::Type{so3})
  w=sqrt(uw'uw)
  return w > eps(Float64) ? so3(uw/w,w) : so3(uw,0)
end

"""
		to_matrix(xr::so3)
"""
function to_matrix(xr::so3)
  (a1,a2,a3)=xr.u*xr.w
  SMatrix{3,3,Float64,9}([0 -a3 a2
                         a3 0 -a1
                         -a2 a1 0])
end 
"""
		to_matrix(Xr::SO3)
"""
to_matrix(Xr::SO3) = Xr.R

# """
#     vee(Rskew::SMatrix{3,3,Float64})
#
# Rskew  -> so3^\\vee = R^3
# """
# vee(Rskew::SMatrix{3,3,Float64}) = [Rskew[3,2],Rskew[1,3],Rskew[2,1]]

"""
		Log(X::SO3) -> SVector3
"""
function Log(X::SO3)
  vee(so3(X.u,X.w))
end

# """
# 		SO3(R::SMatrix{3,3,Float64,9})
#
# # Constructor for ``\\mathrm{SO}(3)``
# # """
# # function SO3(R::SMatrix{3,3,Float64,9})
# #   Rskew=R-R'
# #   @assert is_skew(Rskew)
# #   w=acos((R[1,1]+R[2,2]+R[3,3]-1)/2)
# #   if w > eps(Float64)
# #     uw=0.5*w/sin(w)*SA_F64[Rskew[3,2],Rskew[1,3],Rskew[2,1]]
# #     return SO3(uw/w,w,R)
# #   else
# #     uw=0.5*SA_F64[Rskew[3,2],Rskew[1,3],Rskew[2,1]]
# #     return SO3(uw,w,R)
# #   end
# # end

"""
    Exp(wu::SVector{3,Float64}, ::Union{Type{SO3},Type{so3}})::SO3

``\\mathfrak{so}(3)^\\vee~\\to~\\mathrm{SO}(3)``  
"""
function Exp(wu::SVector{3,Float64}, ::Union{Type{SO3},Type{so3}})
  SO3(wu)
end
"""
    Exp(u::SVector{3,Float64}, w::Float64)::SO3

so3^\\vee -> SO3
"""
function Exp(u::SVector{3,Float64}, w::Float64)
  SO3(u,w)
end
# """
# 		SO3(u::SVector{3,Float64}, w::Float64)
# """
# SO3(u::SVector{3,Float64}, w::Float64) = Exp(u,w)
# """
# 		SO3(uw::SVector{3,Float64})
# """
# SO3(uw::SVector{3,Float64}) = Exp(uw, so3)

"""
		Adjm(Xr::SO3)
"""
function Adjm(Xr::SO3)
  Xr.R
end


"""
    Jr(xr::so3)
"""
function Jr(xr::so3)
  u=xr.u
  w=xr.w
  # defer to the other definition
  if w > eps(Float64)
    return SMatrix{3,3,Float64}([1 0 0;0 1 0;0 0 1]) - (1-cos(w))/w^2*skew(u*w) + (w-sin(w))/w^3*skew(u*w)^2
  else
    # @warn "eps"
    return SMatrix{3,3,Float64}([1 0 0;0 1 0;0 0 1])
  end
end

"""
    Jl(xr::so3)
"""
Jl(xr::so3)=Jr(xr)'

"""
    Jrinv(xr::so3)
"""
function Jrinv(xr::so3)
  u=xr.u
  w=xr.w
  if w > eps(Float64)
    return SMatrix{3,3,Float64}([1 0 0;0 1 0;0 0 1]) + 0.5*skew(u*w) + (1/w^2-(1+cos(w))/(2w*sin(w)))*skew(u*w)^2
  else
    # @warn "eps"
    return SMatrix{3,3,Float64}([1 0 0;0 1 0;0 0 1])
  end
end
"""
    Jlinv(xr::so3)
"""
Jlinv(xr::so3)=Jrinv(xr)'

#--------------------------------------------------------------------#
#                                SE3                                 #
#--------------------------------------------------------------------#
"""
    SE3

Special Euclidean Group 3.

# Fields:
* `t::SVector{3,Float64}`: translation vector
* `rot::SO3`: rotation (Special Orthogonal 3)
"""
struct SE3
  t::SVector{3,Float64}
  rot::SO3

  @doc """
     SE3()
  """
  function SE3()
    new(SA_F64[0,0,0],SO3())
  end
  @doc """
     SE3(t,rot)
  """
  function SE3(t::SVector{3,Float64}, rot::SO3)
    SE3(t,rot)
  end
end

"""
    se3

``\\mathfrak{se}(3)``, Lie algebra of the Special Euclidean Group 3.  Also known as twist.

# Fields
* `v::SVector{3,Float64}`: translational twist.
* `w::so3`: rotational twist.
"""
struct se3
  v::SVector{3,Float64}
  w::so3
  @doc """
      se3(v,uw)
  """
  function se3(v::SVector{3,Float64}, uw::SVector{3,Float64}) 
    new(v, so3(uw))
  end
  @doc """
      se3(v,u,w)
  """
  function se3(v::SVector{3,Float64}, u::SVector{3,Float64},w::Float64)
    new(v, so3(u,w))
  end
end
# se3(v::SVector{3,Float64}, uw::SVector{3,Float64}) = se3(v, so3(uw))
# se3(v::SVector{3,Float64}, u::SVector{3,Float64},w::Float64) = se3(v, so3(u,w))

"""
		vee(x::se3)
"""
vee(x::se3)=SVector{6,Float64}([x.v...,vee(x.w)...])
"""
		hat(x::SVector{6,Float64}) -> se3
"""
hat(x::SVector{6,Float64})=se3(SA_F64[x[1:3]...],SA_F64[x[4:end]...])

"""
		to_matrix(X::SE3)
"""
to_matrix(X::SE3)=SMatrix{4,4,Float64,16}([to_matrix(X.rot) X.t;0 0 0 1])

"""
		Adjm(X::SE3)
"""
function Adjm(X::SE3)
  SMatrix{6,6,Float64,36}(
    [
      X.rot.R   skew(X.t)*X.rot.R
      zeros(3,3) X.rot.R
    ]
  )
end

#  Exp, Log, exp_lie, log_lie, Jl, Jr, Jlinv, Jlinv, QJSE3 (QJSE3 is internal), VTHSO3 (internal)

function internal_VTHSO3(LogX::SVector{3,Float64}) 
  # so3 input
  Jl(uw, so3) # see remark eq (174), Solà 2018
end

"Log(X::SE3) -> SVector6" 
function Log(X::SE3)
  w=Log(X.rot)
  Vinv=Jlinv(hat(w, so3)) # see remark eq (174), Solà 2018
  SA_F64[Vinv*X.t...,w...]
end

"""
		Exp(tau::SVector{6,Float64})::SE3
"""
function Exp(tau::SVector{6,Float64})
  v=SVector{3,Float64}(tau[1:3]...)
  w=SVector{3,Float64}(tau[4:6]...)
  V=Jl(hat(w,so3))  # see remark eq (174), Solà 2018
  SE3(V*v,Exp(w,so3))
end

function internal_QJSE3(vw::se3)
  rho=vw.v
  th=vw.w.w
  wth=vw.w.u*th
  sth=sin(th)
  cth=cos(th)
  thsq=th^2
  thcube=thsq*th
  thquatro=thsq*thsq
  thquint=thquatro*th
  skt=skew(wth)
  skr=skew(rho)
  A=skt*skr
  B=skr*skt
  C=A*skt
  D=skt*A
  E=B*skt
  F=C*skt
  G=skt*C
  if thquint > eps(Float64)
    f3=(th-sth)/thcube
    f4=(1-thsq/2-cth)/thquatro
    f5=(f4-3(th-sth-thcube/6)/thquint)
  else
    # @warn "eps"
    # small angle approximations, based on 
    # manif/impl/se3/SE3Tangent_base.h  l261 commit de52f8b
    f3=1.0/6+thsq/120
    f4=-1.0/24+thsq/720
    f5=-1.0/60
  end
  # TODO: check 0.5 factor in front of f5
  return 0.5skr+f3*(A+B+C)-f4*(D+E-3*C)-0.5*f5*(F+G)
end

"""
		Jl(vw::se3)
"""
function Jl(vw::se3)
  Jlw=Jl(vw.w)
  Q=internal_QJSE3(vw)
  SMatrix{6,6,Float64,36}([Jlw Q; zeros(3,3) Jlw])
end
"""
		Jlinv(vw::se3)
"""
function Jlinv(vw::se3)
  Jlinvw=Jlinv(vw.w)
  Q=internal_QJSE3(vw)
  SMatrix{6,6,Float64,36}([Jlinvw -Jlinvw*Q*Jlinvw; zeros(3,3) Jlinvw])
end

"""
		Jr(vw::se3)
"""
Jr(vw::se3)=Jl(hat(-vee(vw)))
"""
		Jrinv(vw::se3)
"""
Jrinv(vw::se3)=Jlinv(hat(-vee(vw)))

# end SE3/SO3 TODO: move to bottom

# TODO: define \approx
#
# TODO: test (vee(x) |> hat) \approx x
# TODO: test Jlinv(x) \approx inv(Jl(x))

#--------------------------------------------------------------------#
#                                SO2                                 #
#--------------------------------------------------------------------#
"""
		SO2

# Fields
* `th::Float64`: angle (not necessarily inside ``[-\\pi,\\pi]``)
* `c::Float64`:  cos
* `s::Float64`:  sin
"""
struct SO2
    th::Float64 # not necessarily inside [-pi,pi]
    c::Float64 # cos
    s::Float64 # sin

    @doc """
        SO2()
    """
    function SO2()
      new(0,1,0)
    end
    @doc """
        SO2(th::Number)
    """
    function SO2(th::Number)
      sth = sin(th)
      cth = cos(th)
      new(th, cth, sth)
    end
    @doc """
        SO2(co,si)
    """
    function SO2(co,si)
      new(atan(si, co), co, si)
    end
    @doc """
        SO2(R::SMatrix{2,2,Float64,4})
    """
    function SO2(R::SMatrix{2,2,Float64,4})
      new(R[1, 1], R[2, 1]) # TODO: check that R is legit orthogonal
    end
end

"""
    so2

# Fields:
* `w::Float64`: rotational twist
"""
struct so2
    w::Float64

    @doc """
        so2(w::Number=0)
    """
    function so2(w::Number=0)
      new(w)
    end
end

# TODO: CONTINUE HERE for docs

"""
    SE2 { t::SVector{2,Float64}, rot::SO2 }

Possible constructors: SE2(x,y,th), SE2(t, rot), SE2(t, th), SE2(X::SMatrix)
"""
struct SE2
    t::SVector{2,Float64}
    rot::SO2
end
SE2(x, y, th) = begin
    # @info "less efficient default method from SE2 (requires conversions)"
    SE2(convert(Float64, x), convert(Float64, y), convert(Float64, th))
end
SE2(x::Float64, y::Float64, th::Float64) = SE2(SA_F64[x, y], SO2(th))
SE2(t::Vector{Float64}, rot::SO2) = begin
    if length(t) != 2
        throw(DimensionMismatch)
    else
        # @info "SE2 ctor: dynamic vector input"
        SE2(SA_F64[t[1], t[2]], rot)
    end
end
SE2(t::Vector, th) = begin
    @assert(length(t) == 2)
    # @info "less efficient default method from SE2 (requires conversion to vector{Float64})"
    SE2(convert(Vector{Float64}, t), SO2(th))
end
SE2(t::Vector{Float64}, th::Float64) = SE2(t, SO2(th))
SE2(X::SMatrix{3,3,Float64}) = SE2(X[1:2, 3], SO2(X[1:2, 1:2]))
SE2(X::Matrix) = begin
    @assert(size(X) == (3, 3))
    SE2(X[1:2, 3], SO2(X[1:2, 1:2]))
end


"""
    se2 {  vx::Float64, vy::Float64, w::Float64  }
"""
struct se2
    # a.k.a. screw (se3 would be the twist)
    vx::Float64
    vy::Float64
    w::Float64
end
se2(tau::SVector{3,Float64}) = se2(tau...)
se2(tau::Vector{Float64}) = begin
    if length(tau) != 3
        throw(DimensionMismatch)
    else
        # @info "se2 ctor: dynamic vector input"
        se2(tau...)
    end
end
se2(tau) = begin
    @assert length(tau) == 3
    # @info "less efficient default method from se2 (requires conversion to vector{Float64})"
    se2(convert(Vector{Float64}, tau))
end

"""
    ecpi

Angle value of SO2 object, given in bounds [-π,π]
"""
ecpi(r::SO2) = atan(r.s, r.c)
ecpi(a::Float64) = atan(sin(a), cos(a))
ecpi(a) = begin
    # @info "less efficient method for ecpi (requires conversion to float64)"
    ecpi(convert(Float64, a))
end
# force representation between [-pi,pi]; r.th may be outside (use case: external readability)

"""
		to_matrix(r::SO2)
"""
to_matrix(r::SO2) = SA_F64[r.c -r.s; r.s r.c]
"""
		to_matrix(X::SE2)
"""
to_matrix(X::SE2) = begin
    R = to_matrix(X.rot)
    # @info typeof(R)
    SA_F64[
        R[1, 1] R[1, 2] X.t[1]
        R[2, 1] R[2, 2] X.t[2]
        0 0 1
    ]
end
"""
		to_matrix(w::so2)
"""
to_matrix(w::so2) = SA_F64[0 -w.w; w.w 0]
"""
		to_matrix(sk::se2)
"""
to_matrix(sk::se2) = begin
    SA_F64[
        0 -sk.w sk.vx
        sk.w 0 sk.vy
        0 0 0
    ]
end

"""
		vee(w::so2)
"""
vee(w::so2) = w.w
"""
		vee(sk::se2)
"""
vee(sk::se2) = SA_F64[sk.vx, sk.vy, sk.w]

"""
		hat(w::Float64) -> so2
"""
hat(w::Number) = so2(w)
"""
		hat(tau::SVector{3,Float64})::se2
"""
hat(tau::SVector{3,Float64}) = se2(tau...)
hat(tau::Vector{Float64}) = begin
    # @info "hat ctor: dynamic vector input"
    se2(tau...)
end

import Base: *
"""
    *(rot1::SO2, rot2::SO2)
"""
function Base.:*(rot1::SO2, rot2::SO2)
    SO2(rot1.th + rot2.th)
end
"""
    *(rot::SO2, t::SVector{2,Float64})
"""
function Base.:*(rot::SO2, t::SVector{2,Float64}) # action SO2*point
    to_matrix(rot) * t
end
"""
    *(pose::SE2, t::SVector{2,Float64}) 
"""
function Base.:*(pose::SE2, t::SVector{2,Float64}) # action SE2*point
    pose.t + pose.rot * t
end
"""
    *(X1::SE2, X2::SE2)
"""
function Base.:*(X1::SE2, X2::SE2)
    SE2(X1.t + X1.rot * X2.t, X1.rot * X2.rot)
end
"""
    *(Xr1::SO3,Xr2::SO3)
"""
function Base.:*(Xr1::SO3,Xr2::SO3)
  SO3(Xr1.R*Xr2.R)
end
"""
    *(X1::SE3, X2::SE3)
"""
function Base.:*(X1::SE3, X2::SE3)
    SE3(X1.t + X1.rot * X2.t, X1.rot * X2.rot)
end
"""
    *(rot::SO3, t::SVector{3,Float64})
"""
function Base.:*(rot::SO3, t::SVector{3,Float64}) # action SO3*point
    rot.R * t
end
"""
    *(rot1::SO2, rot2::SO2)
"""
function Base.:*(pose::SE3, t::SVector{3,Float64}) # action SE3*point
    pose.t + pose.rot * t
end

import Base: inv
"""
    inv(rot::SO2)
"""
Base.:inv(rot::SO2) = SO2(-rot.th, rot.c, -rot.s)
"""
    inv(X::SE2)
"""
Base.:inv(X::SE2) = SE2(-(inv(X.rot) * X.t), inv(X.rot))
"""
    inv(rot::SO3)
"""
Base.:inv(rot::SO3) = SO3(-rot.u, rot.w, rot.R')
"""
    inv(rot::SE3)
"""
Base.:inv(X::SE3) = SE3(-(inv(X.rot)*X.t),inv(X.rot))  # test: to_matrix(X*inv(X))  ≈ I(4)

"""
    Adjm(rot::SO2)

The adjoint matrix.
"""
Adjm(rot::SO2) = 1

"""
		Adjm(X::SE2)
"""
function Adjm(X::SE2)
    SA_F64[
        X.rot.c -X.rot.s X.t[2]
        X.rot.s X.rot.c -X.t[1]
        0 0 1
    ]
end

"""
		exp_lie(w::so2)
"""
exp_lie(w::so2) = SO2(w.w)

"""
		exp_lie(sk::se2) -> (; SE2, K1, w_sq)
"""
function exp_lie(sk::se2)
    # technique from manifcpp
    # credit: sola/deray & contributors of https://github.com/artivis/manif
    w_sq = sk.w * sk.w
    # K1 is  sin(w)÷w
    # K2 is  (1-cos(w))÷w
    if w_sq < eps(Float64) * 10e9 # WARN_EPS
        K1 = 1 - (1.0 / 6.0) * w_sq
        K2 = 0.5 * sk.w - (1.0 / 24.0) * sk.w * w_sq
    else
        K1 = sin(sk.w) / sk.w
        K2 = (1 - cos(sk.w)) / sk.w
    end
    t = SA_F64[K1 -K2; K2 K1] * SA_F64[sk.vx; sk.vy]
    # returns a named tuple
    (SE2 = SE2(t, Exp(sk.w)), K1 = K1, K2 = K2, w_sq = w_sq)
end

"""
    Exp(w::Float64)::SO2

so2^\\vee -> SO2
"""
Exp(w::Float64) = exp_lie(hat(w))
"""
    Exp(tau::SVector{3,Float64})::SE2

se2^\\vee -> SE2
"""
Exp(tau::SVector{3,Float64}) = exp_lie(hat(tau)).SE2
Exp(tau::Vector{Float64}) = begin
    # @info "Exp function call: dynamic input vector"
    Exp(SA_F64[tau...])
end

"""
    ExpAndJr(tau::SVector{3,Float64}) -> (SE2, Jr_exp)

Returns the group manifold object directly from a vector representation of a Lie algebra element.
Also returns the right Jacobian, as it allows to save some computations compared to doing the 2 things
separately.
"""
function ExpAndJr(tau::SVector{3,Float64})
    # technique from manifcpp
    # credit: sola/deray & contributors of https://github.com/artivis/manif
    expmap_value, K1, K2, w_sq = exp_lie(hat(tau))
    vx = tau[1]
    vy = tau[2]
    w = tau[3]
    if w_sq < eps(Float64) * 10e9 # WARN_EPS
        J13 = -vy / 2.0 + w * vx / 6.0
        J23 = vx / 2.0 + w * vy / 6.0
    else
        sw = sin(w)
        cw = cos(w)
        J13 = (-vy + w * vx + vy * cw - vx * sw) / w_sq
        J23 = (vx + w * vy - vx * cw - vy * sw) / w_sq
    end
    Jr_exp = SA_F64[K1 K2 J13; -K2 K1 J23; 0 0 1]
    expmap_value, Jr_exp
end

"""
		log_lie(rot::SO2) -> so2
"""
log_lie(rot::SO2) = so2(atan(rot.s, rot.c)) # non-bijectivity
"""
		Log(rot::SO2) -> SVector3
"""
Log(rot::SO2) = vee(log_lie(rot))

"""
		log_lie(X::SE2) -> se2
"""
function log_lie(X::SE2)
    # technique from manifcpp
    # credit: sola/deray & contributors of https://github.com/artivis/manif
    th = ecpi(X.rot)
    th_sq = th * th
    # K1 is  sin(w)÷w
    # K2 is  (1-cos(w))÷w
    if th_sq < eps(Float64) * 10e9 # WARN_EPS
        K1 = 1 - (1.0 / 6) * th_sq
        K2 = 0.5 * th - (1.0 / 24.0) * th * th_sq
    else
        K1 = sin(th) / th
        K2 = (1 - cos(th)) / th
    end
    Vinv = SA_F64[K1 K2; -K2 K1] / (K1 * K1 + K2 * K2)
    se2(Vinv * X.t..., th)
end

"""
		Log(X::SE2) -> SVector3
"""
Log(X::SE2) = vee(log_lie(X))

"""
		Jrinv
"""
function Jrinv(sk::se2)
    # technique from manifcpp
    # credit: sola/deray & contributors of https://github.com/artivis/manif
    vx = sk.vx
    vy = sk.vy
    w = sk.w
    sw = sin(w)
    cw = cos(w)
    wsw = w * sw
    wcw = w * cw
    w_sq = w * w
    J12 = -w / 2
    J21 = -J12
    if (w_sq > eps(Float64) * 10e9) # WARN_EPS
        J11 = -wsw / (2 * cw - 2)
        J22 = J11
        dd = 2 * w * (cw - 1)
        J13 = (wsw * vx + wcw * vy - w * vy + 2 * vx * cw - 2 * vx) / dd
        J23 = (-wcw * vx + wsw * vy + w * vx + 2 * vy * cw - 2 * vy) / dd
    else
        J11 = 1.0 - w_sq / 12.0
        J22 = J11
        J13 = vy / 2 + w * vx / 12
        J23 = -vx / 2 + w * vy / 12
    end
    J31 = 0
    J32 = 0
    J33 = 1
    SA_F64[J11 J12 J13; J21 J22 J23; J31 J32 J33]
end
function Jrinv(isomorph_to_se2::SVector{3,Float64})
    Jrinv(se2(isomorph_to_se2...))
end

"""
		Jr
"""
function Jr(sk::se2)
    # technique from manifcpp
    # credit: sola/deray & contributors of https://github.com/artivis/manif
    # K1 is  sin(w)÷w
    # K2 is  (1-cos(w))÷w
    vx = sk.vx
    vy = sk.vy
    w = sk.w
    w_sq = w * w
    if w_sq < eps(Float64) * 10e9 # WARN_EPS
        K1 = 1 - w_sq / 6
        K2 = 0.5 * sk.w - 1.0 / 24 * sk.w * w_sq
        J13 = -vy / 2.0 + w * vx / 6.0
        J23 = vx / 2.0 + w * vy / 6.0
    else
        sw = sin(w)
        cw = cos(w)
        K1 = sw / sk.w
        K2 = (1 - cw) / sk.w
        J13 = (-vy + w * vx + vy * cw - vx * sw) / w_sq
        J23 = (vx + w * vy - vx * cw - vy * sw) / w_sq
    end
    SA_F64[K1 K2 J13; -K2 K1 J23; 0 0 1]
end

"""
		Jlinv
"""
Jlinv(sk::se2) = Jrinv(hat(-vee(sk)))

"""
		Jl
"""
Jl(sk::se2) = Jr(hat(-vee(sk)))

"""
    wedge

Alias for hat.
"""
wedge=hat

# TODO: to_quat for SO3, to_S1 (unit circle) for SO2
# TODO: approx check between several SE2 se2 SO2 so2

end

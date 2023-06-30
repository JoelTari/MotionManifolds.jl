module ManifoldExtras

using StaticArrays

export SO2, so2, SE2, se2, to_matrix, hat, vee, Adjm, ecpi, Log, Exp, log_lie, exp_lie, Jr, Jrinv, Jl, Jlinv, ExpAndJr

export SO2FromMat

"SO2"
struct SO2
  th::Float64 # not necessarily inside [-pi,pi]
  c::Float64 # cos
  s::Float64 # sin
end
# generic boilerplate stuff,
# but you get lectured with logs
SO2(th)=begin
  @info "less efficient default method, non float64 scalar input"
  sth = sin(th)
  cth = cos(th)
  SO2(th,cth,sth)
end 
SO2(co,si)=begin
  @info "less effficient default method for si,co ctor"
  SO2(atan(si, co), co, si)
end
# prefer methods
SO2(th::Float64) = begin
  sth = sin(th)
  cth = cos(th)
  SO2(th,cth,sth)
end
SO2(cosinus::Float64, sinus::Float64) = SO2(atan(sinus, cosinus), cosinus, sinus)
SO2(R::SMatrix{2,2,Float64,4}) = SO2(R[1,1],R[2,1]) # TODO: check that R is legit orthogonal
SO2FromMat(R::Matrix{Float64}) = begin
  if size(R) != (2,2)
    throw(DimensionMismatch)
  else
    @info "SO2 ctor: dynamic matrix input"
    SO2(atan(R[1,1],R[2,1]))
  end
end
SO2FromMat(R::Matrix)=begin
  @info "less efficient default method for R matrix"
  if size(R) != (2,2)
    throw(DimensionMismatch)
  else
    @info "SO2 ctor: dynamic matrix input"
    SO2(atan(R[1,1],R[2,1]))
  end
end

"SE2"
struct SE2
  t::SVector{2,Float64}
  rot::SO2
end
SE2(x,y,th)= begin
  @info "less efficient default method from SE2 (requires conversions)"
  SE2(convert(Float64,x),convert(Float64,y),convert(Float64,th))
end
SE2(x::Float64,y::Float64,th::Float64) = SE2(SA_F64[x,y],SO2(th))
SE2(t::Vector{Float64},rot::SO2) = begin
  if length(t) != 2
    throw(DimensionMismatch)
  else
    @info "SE2 ctor: dynamic vector input"
    SE2(SA_F64[t[1],t[2]],rot)
  end
end
SE2(t::Vector, th)=begin
  @assert(length(t)==2)
  @info "less efficient default method from SE2 (requires conversion to vector{Float64})"
  SE2(convert(Vector{Float64}, t), SO2(th))
end
SE2(t::Vector{Float64},th::Float64) = SE2(t,SO2(th))
SE2(X::SMatrix{3,3,Float64}) = SE2(X[1:2,3], SO2(X[1:2,1:2]))
SE2(X::Matrix) = begin
  @assert(size(X)==(3,3))
  SE2(X[1:2,3], SO2(X[1:2,1:2]))
end

"so2"
struct so2
  w::Float64
end

"se2"
struct se2
  # a.k.a. screw (se3 would be the twist)
  vx::Float64
  vy::Float64
  w::Float64
end
se2(tau::SVector{3,Float64}) = se2(tau...)
se2(tau::Vector{Float64}) = begin
  if length(tau)!=3
    throw(DimensionMismatch)
  else
    @info "se2 ctor: dynamic vector input"
    se2(tau...)
  end
end
se2(tau)=begin
  @assert length(tau)==3
  @info "less efficient default method from se2 (requires conversion to vector{Float64})"
  se2(convert(Vector{Float64},tau))
end

"""
    ecpi

Angle value of SO2 object, given in bounds [-π,π]
"""
ecpi(r::SO2) = atan(r.s,r.c) 
ecpi(a::Float64) = atan(sin(a),cos(a)) 
ecpi(a) = begin
  @info "less efficient method for ecpi (requires conversion to float64)"
  ecpi(convert(Float64,a))
end
# force representation between [-pi,pi]; r.th may be outside (use case: external readability)

"to_matrix"
to_matrix(r::SO2) = SA_F64[r.c -r.s;r.s r.c]
to_matrix(X::SE2) = begin
  R = to_matrix(X.rot)
  @info typeof(R)
  SA_F64[R[1,1] R[1,2] X.t[1]
         R[2,1] R[2,2] X.t[2]
         0      0      1    ]
end
to_matrix(w::so2) = SA_F64[0 -w.w;w.w 0]
to_matrix(sk::se2) = begin
  SA_F64[0       -sk.w    sk.vx
         sk.w    0        sk.vy
         0       0        0    ]
end

"vee"
vee(w::so2) = w.w
vee(sk::se2) = SA_F64[sk.vx,sk.vy,sk.w]

"hat"
hat(w::Float64) = so2(w)
hat(tau::SVector{3,Float64}) = se2(tau...)
hat(tau::Vector{Float64}) = begin
  @info "hat ctor: dynamic vector input"
  se2(tau...)
end

import Base: *
"*"
function Base.:*(rot1::SO2,rot2::SO2)  
  SO2(rot1.th+rot2.th)
end
function Base.:*(rot::SO2, t::SVector{2,Float64}) # action SO2*point
  to_matrix(rot)*t
end
function Base.:*(X1::SE2,X2::SE2)
  SE2(X1.t + X1.rot*X2.t,X1.rot*X2.rot)
end

import Base: inv
"inv"
Base.:inv(rot::SO2) = SO2(-rot.th,rot.c,-rot.s)
Base.:inv(X::SE2) = SE2(-(inv(X.rot)*X.t),inv(X.rot))

"""
    Adjm

The adjoint matrix.
"""
Adjm(rot::SO2) = 1

function Adjm(X::SE2)
SA_F64[X.rot.c -X.rot.s    X.t[2]
       X.rot.s  X.rot.c   -X.t[1]
       0        0          1      ]
end

"exp_lie"
exp_lie(w::so2)=SO2(w.w)

# """
#     exp_lie
#
# Returns the Group Manifold object from a Lie algebra element.
# """
function exp_lie(sk::se2)
  # technique from manifcpp
  # credit: sola/deray & contributors of https://github.com/artivis/manif
  w_sq = sk.w*sk.w    
  # K1 is  sin(w)÷w 
  # K2 is  (1-cos(w))÷w     
  if w_sq < eps(Float64)*10e9 # WARN_EPS
    K1 = 1 - (1.0/6.0)*w_sq
    K2 = .5*sk.w - (1.0/24.0)*sk.w*w_sq
  else
    K1 = sin(sk.w)/sk.w;
    K2 = (1-cos(sk.w))/sk.w;
  end
  t = SA_F64[K1 -K2;K2 K1]*SA_F64[sk.vx;sk.vy];
  SE2(t,Exp(sk.w)), K1, K2, w_sq
end

"Exp"
Exp(w::Float64) = exp_lie(hat(w))
Exp(tau::SVector{3,Float64}) = exp_lie(hat(tau))[1]
Exp(tau::Vector{Float64}) = begin
  @info "Exp function call: dynamic input vector"
  Exp(SA_F64[tau...])
end

"""
    ExpAndJr

Returns the group manifold object directly from a vector representation of a Lie algebra element.
Also returns the right Jacobian, as it allows to save some computations compared to doing the 2 things
separately.
"""
function ExpAndJr(tau::SVector{3,Float64})
  # technique from manifcpp
  # credit: sola/deray & contributors of https://github.com/artivis/manif
  expmap_value, K1 , K2, w_sq = exp_lie(hat(tau))
  vx = tau[1]; vy=tau[2]; w=tau[3]
  if w_sq < eps(Float64)*10e9 # WARN_EPS
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

"log_lie"
log_lie(rot::SO2) = so2(atan(rot.s,rot.c)) # non-bijectivity
Log(rot::SO2)=vee(log_lie(rot))

function log_lie(X::SE2)
  # technique from manifcpp
  # credit: sola/deray & contributors of https://github.com/artivis/manif
  th=ecpi(X.rot)
  th_sq = th*th
  # K1 is  sin(w)÷w 
  # K2 is  (1-cos(w))÷w     
  if th_sq < eps(Float64)*10e9 # WARN_EPS
    K1 = 1 - (1.0/6)*th_sq
    K2 = .5*th - (1.0/24.0)*th*th_sq 
  else
    K1 = sin(th)/th;
    K2 = (1-cos(th))/th;
  end
  Vinv =SA_F64[K1 K2;-K2 K1]/(K1*K1+K2*K2)
  se2(Vinv*X.t ..., th)
end

"Log"
Log(X::SE2) = vee(log_lie(X))

"Jrinv"
function Jrinv(sk::se2)
  # technique from manifcpp
  # credit: sola/deray & contributors of https://github.com/artivis/manif
  vx = sk.vx; vy=sk.vy; w=sk.w
  sw = sin(w)
  cw = cos(w)
  wsw = w*sw
  wcw = w*cw
  w_sq = w*w
  J12 = -w/2
  J21 = -J12
  if (w_sq > eps(Float64)*10e9) # WARN_EPS
    J11= - wsw/(2*cw-2)
    J22=J11
    dd=2*w*(cw-1)
    J13=(wsw*vx+ wcw*vy -w*vy + 2*vx*cw -2*vx)/dd
    J23=(-wcw*vx+wsw*vy+w*vx+ 2*vy*cw-2*vy)/dd
  else
    J11=1.0-w_sq/12.0
    J22=J11
    J13 = vy/2+ w*vx/12
    J23 = -vx/2+ w*vy/12
  end
  J31=0
  J32=0
  J33=1
  SA_F64[J11 J12 J13;J21 J22 J23;J31 J32 J33]
end

"Jr"
function Jr(sk::se2)
  # technique from manifcpp
  # credit: sola/deray & contributors of https://github.com/artivis/manif
  # K1 is  sin(w)÷w 
  # K2 is  (1-cos(w))÷w     
  vx = sk.vx; vy=sk.vy; w=sk.w
  w_sq = w*w
  if w_sq < eps(Float64)*10e9 # WARN_EPS
    K1 = 1 - w_sq/6
    K2 = .5*sk.w - 1.0/24*sk.w*w_sq
    J13 = -vy/2.0+w*vx/6.0
    J23 =  vx/2.0+w*vy/6.0
  else
    sw = sin(w)
    cw = cos(w)
    K1 = sw/sk.w;
    K2 = (1-cw)/sk.w;
    J13 = (-vy+w*vx+vy*cw-vx*sw)/w_sq
    J23 = ( vx+w*vy-vx*cw-vy*sw)/w_sq
  end
  SA_F64[K1 K2 J13;-K2 K1 J23;0 0 1]
end

"Jlinv"
Jlinv(sk::se2) = Jrinv( hat(-vee(sk)) )

"Jl"
Jl(sk::se2) = Jr( hat(-vee(sk)) )


# TODO: to_quat for SO3, to_S1 (unit circle) for SO2
# TODO: SO3 SE3 so3 se3
# TODO: approx check between several SE2 se2 SO2 so2

end

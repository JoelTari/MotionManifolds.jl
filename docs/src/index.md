# MotionManifolds.jl

```@meta
DocTestSetup = :(using MotionManifolds, StaticArrays)
```

# Overview

Documentation for MotionManifolds.jl, a small library to manipulate matrix Lie groups
frequently encountered in robotics.

A working minimal example:

```jldoctest
julia> skew(SA_F64[1,2,3])
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.0  -3.0   2.0
  3.0   0.0  -1.0
 -2.0   1.0   0.0
```

Contents:

```@contents
```

## Library

The convention is that Lie groups are denoted with upper case letters, while 
their respective Lie algebra are denoted with lower case letters.

### Special Orthogonal 2 group, a.k.a. SO2

```@docs
SO2
SO2()
SO2(th::Number)
SO2(co,si)
SO2(R::SMatrix{2,2,Float64,4})
```

```@docs
so2
so2(w::Number=0)
```
### Special Euclidean 2 group, a.k.a. SE2

```@docs
SE2
SE2()
SE2(t::Vector{Float64}, rot::SO2)
SE2(x::Number, y::Number, th::Number)
SE2(X::SMatrix{3,3,Float64})
```

```@docs
se2
se2(vx,vy,w)
se2()
se2(tau::SVector{3,Float64})
```
### Special Orthogonal 3 group, a.k.a. SO3

```@docs
SO3
SO3()
SO3(u::SVector{3,Float64},w::Float64,R::SMatrix{3,3,Float64,9})
SO3(R::SMatrix{3,3,Float64,9})
SO3(u::SVector{3,Float64},w::Float64)
SO3(uw::SVector{3,Float64})
```

```@docs
so3
so3()
so3(uw::SVector{3,Float64})
```

### Special Euclidean 3 group, a.k.a. SE3

```@docs
SE3
SE3()
SE3(t::SVector{3,Float64}, rot::SO3)
```

```@docs
se3
se3(v::SVector{3,Float64}, uw::SVector{3,Float64})
se3(v::SVector{3,Float64}, u::SVector{3,Float64},w::Float64)
```


## Manifold manipulations

```@docs
*
inv
```

```@docs
hat
wedge
vee
```

```@docs
exp_lie
log_lie
Exp
Log
Adjm
```

## Computation of derivatives

```@docs
Jr
Jrinv
Jl
Jlinv
```

```@docs
ExpAndJr
```

## Utils

```@docs
skew
is_skew
to_matrix
ecpi
```

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

### 

Special Orthogonal 2:
```@docs
SO2
so2
```
Special Euclidian 2
```@docs
SE2
se2
```
### Special Orthogonal 3

Special Orthogonal 3:
```@docs
SO3
so3
so3(uw::SVector{3,Float64})
```
Special Euclidian 3
```@docs
SE3
se3
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

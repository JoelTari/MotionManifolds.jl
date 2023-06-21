module ManifoldExtras

using Manifolds
using RecursiveArrayTools
using StaticArrays
# add another dependency to confirm the workflow
# RecursiveArrayTools
# StaticArrays

export greet, SE2, SE3, SO2, SO3, NV, makeSO2, is_point, angleSO2, makeSE2, getTranslationComponent, getRotationComponent, SE2LieGroup, SE3LieGroup, SO2LieGroup, SO3LieGroup, angleSE2

"""
    greet()

Return a string that greets people.

# Example
```julia
julia> greet()
"Hello ManifoldExtras."
```

# To find out more
See also [documentation](https://docs.julialang.org/en/v1/manual/documentation/).
"""
function greet()
  "Hello ManifoldExtras."
end

"The Special Euclidian 2 Manifold: \$\\mathrm{SE}(n)\$"
const SE2LieGroup=SpecialEuclidean(2)
const SE3LieGroup=SpecialEuclidean(2)
const SO2LieGroup=SpecialOrthogonal(2)
const SO3LieGroup=SpecialOrthogonal(2)

"""
    NV

NV is the number of nodes ``n``, valued 6.

# Example
```jldoctest; setup = :(using ManifoldExtras)
julia> NV
6
```
"""
const NV = 6

const is_point = Manifolds.is_point # TODO: export Manifolds instead

"Type alias"
const SE2=ArrayPartition{Float64, Tuple{SVector{2, Float64}, SMatrix{2,2,Float64, 4}}}
const SE3=ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}}
const SO2=SMatrix{2,2,Float64,4}
const SO3=SMatrix{3,3,Float64,9}

"Construct a SO2 object, compatible with Manifolds.jl"
function makeSO2(a::Number)
  aa=Float64(a)
  SA_F64[cos(aa) -sin(aa);sin(aa) cos(aa)]
end

# TODO: makeSO3, makeSE3

"get angle from an SO2 object, compatible with Manifolds.jl"
function angleSO2(R::Union{Matrix{Float64},SMatrix{2,2,Float64}})
  atan(R[2,1],R[1,1])
end

"extract translation vector from SE2"
function getTranslationComponent(X::SE2)
  SA_F64[X[1:2]...]
end

"extract rotation matrix from SE2"
function getRotationComponent(X::SE2)
  # reshape(SA_F64[X[3:end]...],2,2)
  SMatrix{2,2,Float64,4}([X[3:end]...])
end

# "extract translation vector from SE3"
# function getTranslationComponent(X::ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}})
# end
# "extract rotation matrix from SE3"
# function getRotationComponent(X::ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3,3,Float64, 9}}})
#   reshape(SA_F64[X[3:end]...],2,2)
# end


"get angle from an SE2 object, compatible with Manifolds.jl"
function angleSE2(X::SE2)
  angleSO2(getRotationComponent(X))
end


"Construct an SE2 object from x,y,a."
function makeSE2(x::Number,y::Number,a::Number)
  ArrayPartition(SA_F64[Float64(x),Float64(y)],makeSO2(Float64(a)))
end

"Construct an SE2 object from (t, angle), compatible with Manifolds.jl"
function makeSE2(t,a::Number)
  ArrayPartition(SA_F64[t...],makeSO2(Float64(a)))
end

"Construct an SE2 object from (t, R), compatible with Manifolds.jl"
function makeSE2(t,R::SO2)
  ArrayPartition(SA_F64[t...],R)
end

# function 
#
# end
#
# hat
# vee
# Exp
# Log
# Adjm adjoint matrix
#
# rpySE3 quatSE3

end

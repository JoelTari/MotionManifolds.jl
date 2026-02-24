# MotionManifolds

`MotionManifolds` is a Julia package to compute Special Orthogonal and/or Special Euclidean objects, i.e., $\text{SO}(2)$, $\text{SO}(3)$ $\text{SE}(2)$ and $\text{SE}(3)$.

## Quickstart

```julia
git clone https://github.com/JoelTari/MotionManifolds.jl MotionManifolds
cd MotionManifolds
julia --project=. -e "using Pkg; Pkg.precompile()"
```

## Usage

Declare objects:
```julia
julia> using MotionManifolds

julia> x = 1.0; y = 0.0; theta = pi/3;

julia> X = SE2(x, y, theta) # declare SE2 object
SE2([1.0, 0.0], SO2(1.0471975511965976, 0.5000000000000001, 0.8660254037844386))

julia> to_matrix(X)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 0.5       -0.866025  1.0
 0.866025   0.5       0.0
 0.0        0.0       1.0

julia> Y = rand(SO3); # random SO3 object

julia> to_quat(Y)
Quaternion(0.7864336220880247, 0.44836802694679684, 0.3699299570067636, 0.2089021239009323)

julia> y = Log(Y) # projection to 
3-element SVector{3, Float64} with indices SOneTo(3):
 0.966578657184116
 0.7974841639147566
 0.45034507874718005

julia> yhat = hat(y, so3) # so3 object

julia> vee(yhat) ≈  y
true
```

Derivatives and ajoints:
```julia
julia> x = rand(se3);

julia> Jr(x)
6×6 SMatrix{6, 6, Float64, 36} with indices SOneTo(6)×SOneTo(6):
  0.944793     0.0664991  -0.27411   -0.112466     0.33391    -0.202541
 -0.00637613   0.982418    0.161118  -0.223876    -0.0452179   0.195762
  0.281988    -0.146893    0.929075   0.270563    -0.0732738  -0.12895
  0.0          0.0         0.0        0.944793     0.0664991  -0.27411
  0.0          0.0         0.0       -0.00637613   0.982418    0.161118
  0.0          0.0         0.0        0.281988    -0.146893    0.929075

julia> Jl(x) ≈ Jr(-x)
true

julia> Adjm(Exp(x))
6×6 SMatrix{6, 6, Float64, 36} with indices SOneTo(6)×SOneTo(6):
  0.836816  0.0186794   0.547166  -0.327462  -0.37169     0.513498
  0.159037  0.948028   -0.27559    0.691523  -0.132075   -0.0552746
 -0.523877  0.317637    0.790354  -0.313142   0.416052   -0.374771
  0.0       0.0         0.0        0.836816   0.0186794   0.547166
  0.0       0.0         0.0        0.159037   0.948028   -0.27559
  0.0       0.0         0.0       -0.523877   0.317637    0.790354

julia> Adjm(Exp(x)) ≈ Jl(x)*Jrinv(x)  # \mathbf{Ad}_{\text{Exp}(x)} = \mathbf{J}_l(x) \mathbf{J}_r^{-1}(x)
true
```

Operations:
```julia
julia> X = rand(SE3); x = rand(se3);

julia> X + x  # X \circ Exp(x)
SE3([1.3783088969943498, 1.8731290028869272, 1.220467938988993], SO3([0.589074838514794, 0.6558178263889646, 0.47211631323141456], 1.3382375758717824, [0.49750275127048393 -0.16211688643810046 0.8521790466855631; 0.7566968754870624 0.5614415033141072 -0.3349526488692279; -0.4241472045250518 0.8114810863402863 0.40199203400760647]))

julia> Y = rand(SE3);

julia> X - Y == X + inv(Y)  #  
true

julia> b = rand(SVector{3, Float64})

julia> X + b    # action   X.b
```

## Documentation

Documentation of the API. 

## Licence

BSD 2.0

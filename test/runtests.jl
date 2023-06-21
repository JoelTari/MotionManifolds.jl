using ManifoldExtras
using Test

@testset "ManifoldExtras.jl" begin
    # Write your tests here.
    x=-2.3;y=6.5;a=-5*pi/6
    X = makeSE2(x,y,a)
    @test ManifoldExtras.greet() == "Hello ManifoldExtras."
    @test NV==6
    @test [cos(a) -sin(a); sin(a) cos(a)] ≈ getRotationComponent(X)        
end

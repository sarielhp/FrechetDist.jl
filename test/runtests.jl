using FrechetDist
using FrechetDist.cg
using Test

function  example_01()
    P = Polygon2F( );
    Q = Polygon2F( );

    push!( P, point( 0.0, 0 ), point( 1.0, 0.0 ) );
    push!( Q, point( 0.0, 1.0 ), point( 1.0, 1.0 ) );
    return  P, Q;
end

@testset "FrechetDist.jl" begin

    # Write your tests here.
    P,Q = example_01();

    @test FrechetDist.frechet_c_compute( P, Q ).leash == 1.0;
end

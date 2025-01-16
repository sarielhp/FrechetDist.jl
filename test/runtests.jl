using FrechetDist
using FrechetDist.cg: Polygon2F, Point2F
using FrechetDist.cg.polygon
using Test

function  example_01()
    P = Polygon2F( );
    Q = Polygon2F( );

    push!( P, Point2F( 0.0, 0 ), Point2F( 1.0, 0.0 ) );
    push!( Q, Point2F( 0.0, 1.0 ), Point2F( 1.0, 1.0 ) );
    return  P, Q;
end


function  example_02()
    P = Polygon_random( 3, Float64, 100 );
    Q = Polygon_random( 3, Float64, 80 );

    return  P, Q;
end

@testset "FrechetDist.jl" begin

    # Write your tests here.
    P,Q = example_01();

    @test FrechetDist.frechet_c_compute( P, Q ).leash == 1.0;

    P,Q = example_02();

    # Code generation - forcing to generate all different distance
    # calculations
    m_d = frechet_d_compute( P, Q )
    @test m_d.leash > 0.0

    m_d_r = frechet_d_r_compute( P, Q )
    @test m_d_r.leash > 0.0
    @test m_d.leash == m_d_r.leash;

    m_ve = frechet_ve_r_compute( P, Q )
    @test m_ve.leash > 0.0

    m_e = frechet_c_compute( P, Q )
    @test m_e.leash > 0.0

    @test m_e.leash >= m_ve.leash;

    m_a_2 = frechet_c_approx( P, Q, 2.0 );

    @test ( m_a_2.leash > 0.0 );

    @test ( (m_e.leash <= m_a_2.leash) &&  ( m_a_2.leash <= 2.0*m_e.leash ) );

    d = frechet_dist_upper_bound( P, Q );
    @test m_e.leash <= d;

    m_mono, f, P_r, Q_r = frechet_mono_via_refinement( P, Q, 1.01 );
    @test m_mono.leash >= m_ve.leash

    w_P = frechet_width_approx( P );

    w_Q = frechet_width_approx( Q );
    @test   w_P >= 0;
    @test   w_Q >= 0;
end

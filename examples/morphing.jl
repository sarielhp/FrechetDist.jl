#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon
using FrechetDist.cg.segment
using Profile
using InteractiveUtils
using Plots; 

include( "graphics.jl" )

######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T}
     )::Int64  where {D,T}
    ratio::Float64 = 5.0;

    #P = polygon.reverse( P );
    
    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    m_e = frechet_c_compute( P, Q );
    m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    #m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    ma = frechet_c_approx( P, Q, 1.1 );
    println( "Frechet_c_compute:" );
    @time m_e = frechet_c_compute( P, Q );
    println( "Frechet_move_via_refinement:" );
    @time m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    println( "Frechet_c_approx:" );
    @time ma = frechet_c_approx( P, Q, 4.0 );

    n = 400;
    X,Y = morphing_profile( m, n );
    XA,YA = morphing_profile( ma, n );
    XE,YE = morphing_profile( m_e, n );
    
    plt = plot(X, [Y, YA, YE],  dpi=400)
    savefig( plt, "P.png" );
    
    output_frechet_movie_mp4( m, "matching.mp4" );
    
    return  0;
end


function  test_files( fl_a, fl_b )
    println( "---------------------------------------------" );
    println( fl_a );
    println( fl_b );
    println( "---------------------------------------------" );
    P = read_file( fl_a );
    Q = read_file( fl_b );

    frechet_comp( P, Q );

    println( "Done..." );
end



num_args = length( ARGS );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

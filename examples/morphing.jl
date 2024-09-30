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

    m_c = frechet_c_compute( P, Q );
    m_a = frechet_c_approx( P, Q, 1.1 );
    #m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];

    cardi = length( P ) + length( Q ) 
    P_d, Q_d = P, Q;
    if  ( cardi < 100 )
        P_d = Polygon_sample_uniformly( P, 60*length( P ) );
        Q_d = Polygon_sample_uniformly( Q, 60*length( Q ) );
    end
    m_d = frechet_d_compute( P_d, Q_d );


    
    n = 400;
    X_C,Y_C = morphing_profile( m_c, n );
    X_A,Y_A = morphing_profile( m_a, n );
#    XE,YE = morphing_profile( m_e, n );
    X_D,Y_D = morphing_profile( m_d, n );

    plt = plot(X_C, [ Y_C, Y_D],  dpi=400, lw=4,
              labels=["Retractable" "Fréchet"] );
    savefig( plt, "P.png" );

    #output_frechet_movie_mp4( m, "matching.mp4" );

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

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

    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    lP = total_length( P );
    m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    
    PA, QA = Morphing_as_polygons( m );
    output_polygons_to_file( [ PA, QA ], "matching.pdf", false, false, );

    PB,QB,times = Morphing_as_function_w_times( m );
    X = range(0, lP, length=400)
    Y = Vector{Float64}();
    
    for x ∈ X
        p, q = polygons_get_loc_at_time( PB, QB, times, x / lP );
        push!( Y, Dist(p,q ) );
    end

    plt = plot(X, [Y],  dpi=400)
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

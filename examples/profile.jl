#! julia

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

function do_profiles( polys::Vector{Polygon{D,T}} )  where {D,T}
    ratio::Float64 = 5.0;

    M = Vector{Morphing{D,T}}()
    for  i in 2:length( polys )
        push!( M, frechet_c_compute( polys[ 1 ], polys[ i ] ) );
    end

    profiles = Vector{Vector{Float64}}();
    
    n = 400;
    for  i in 1:length( M );
        X,Y = morphing_profile( M[ i ], n );
        ( i == 1 )  &&  push!( profiles, X );
        push!( profiles, Y );
    end

    plt = plot( profiles[1], profiles[2:end],  dpi=400, lw=4,
              labels=["Retractable" "Fréchet"] );
    savefig( plt, "P.png" );

    #output_frechet_movie_mp4( m, "matching.mp4" );

    return  0;
end



num_args = length( ARGS );

polys = Vector{Polygon2F}()

if   num_args > 1
    for  i  in 1:num_args
        println( i );
        P = read_file( ARGS[ i ] );
        ( length( P ) < 10 )  &&  continue;
        
        push!( polys, P );
    end
    do_profiles( polys );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

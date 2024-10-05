#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon
using FrechetDist.cg.segment
using Profile
using InteractiveUtils
using Plots, Printf;

include( "graphics.jl" )

######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T},
                       approx::Float64, prefix 
     )::Int64  where {D,T}
    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    m_a = frechet_c_approx( P, Q, approx );
    println( "F dist in range: ", m_a.lower_bound, " ... ", m_a.leash );
    println( "Ratio: ", m_a.ratio );

    PS, QS =   simplify_morphing_sensitive( m_a, 1.0 );

    println( length( PS ), " = ", length( P ) );
    output_polygons_to_file( [P, PS ], prefix*"curves_P.pdf", true );
    output_polygons_to_file( [Q, QS ], prefix*"curves_Q.pdf", true );
    output_polygons_to_file( [P, Q, PS, QS ], prefix*"together.pdf", true );

    output_polygons_to_file( [P, Q ], prefix*"original.png", false );
    output_polygons_to_file( [P, PS ], prefix*"curves_P.png", false );
    output_polygons_to_file( [Q, QS ], prefix*"curves_Q.png", false );
    output_polygons_to_file( [P, Q, PS, QS ], prefix*"together.png", false );

    n = 100;
    X_A,Y_A = morphing_profile( m_a, n );

    str = @sprintf( "%0.2f", m_a.ratio );
    plt = plot(X_A, [ Y_A ],  dpi=800, lw=4,
              labels="Profile", title=str * " Approximation" );
    #println( Y_A );
    println( maximum( Y_A ) );
    hline!(plt, [m_a.lower_bound], linestyle=:dash, labels="Lower bound")
    savefig( plt, prefix * "P_profile.png" );

    return 0;
end


function  test_files( fl_a, fl_b )
    println( "---------------------------------------------" );
    println( fl_a );
    println( fl_b );
    println( "---------------------------------------------" );
    P = read_file( fl_a );
    Q = read_file( fl_b );

    frechet_comp( P, Q, 4.0, "A_" );
    frechet_comp( P, Q, 1.1, "B_" );
    frechet_comp( P, Q, 1.01, "C_" );

    println( "Done..." );
end



num_args = length( ARGS );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using Profile
using InteractiveUtils

######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T}
     )::Int64  where {D,T}
    ratio::Float64 = 5.0;

    #println( P );
    #println( Q );
    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    #=
    m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    m_d = frechet_c_compute( P, Q );
    =#
    l_min_a, l_max_a =  FEVER_compute_range( P, Q, 10000000.0 );
    #=println( l_min_a, "...", l_min_b );
    @time m = frechet_mono_via_refinement( P, Q, 1.00001 )[1];
    @time m_d = frechet_c_compute( P, Q );
    @time m_d = frechet_c_compute( P, Q );

    println( "frechet_compute range..." );
    @time frechet_ve_r_compute_range( P, Q, 10000000.0 )

    println( "FEVER_compute range..." );
    =#
    @time l_min_b, l_max_b = FEVER_compute_range( P, Q, 10000000.0 );

    println( l_min_a, " ", l_min_b );
    println( l_max_a, " ", l_max_b );
    
    #=
    println( "FEVER_compute range..." );
    @time FEVER_compute_range( P, Q, 10000000.0 );
    println( "FEVER_compute range..." );
    @time FEVER_compute_range( P, Q, 10000000.0 );
    println( m.leash );
    println( m_d.leash );
    =#
    exit( -1 );
    
    for  i in 1:10
        println( "approx( #", cardin(P ), ", ", cardin(Q), ")   approx: ",
                 ratio );
        m = frechet_c_approx( P, Q, ratio );
        exit( -1 );
        if  ( m.ratio == 1.0 )
            return  0;
        end
        lb = m.leash / m.ratio;

        #println( "m.leash: ", m.leash );
        ratio = (ratio - 1.0) / 6.0 + 1.0; # min( m.ratio, 1.01 );
        ratio = min( ratio, 1.1 );
#        println( "ratio: ", ratio );
        if  ( ratio <= 11.0 )
            println( "frechet_c_compute..." );
            m = frechet_c_compute( P, Q );
#            println( "... done" );
            println( "----------------------------------------------" );
            println( "DISTANCE: ", m.leash );
            return  0;
        end
#        if  i > 2
#            println( "RATIO  ", i, " : ", ratio );
#        end
    end
    println( "UNDECIDED" );

    return  0;
end


function  test_files( fl_a, fl_b )
    println( "---------------------------------------------" );
    println( fl_a );
    println( fl_b );
    println( "---------------------------------------------" );
    P = Polygon_read_file( fl_a );
    Q = Polygon_read_file( fl_b );

    #@profile
    #@code_warntype
    #print( @report_opt frechet_comp( P, Q ) );
    #frechet_ve_r_compute_mono_dist( P, Q );

    frechet_comp( P, Q );
    #Profile.print();
    
    println( "Done..." );
    #println( "Sign: ", sgn, "\n\n\n\n\n\n" );
end



num_args = length( ARGS );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

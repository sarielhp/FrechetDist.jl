#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using Profile
using InteractiveUtils
using JET

######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T}
     )::Int64  where {D,T}
    ratio::Float64 = 5.0;

    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    for  i in 1:10
        println( "approx( #", cardin(P ), ", ", cardin(Q), ")   approx: ",
                 ratio );
        m = frechet_c_approx( P, Q, ratio );
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
    print( @report_opt frechet_comp( P, Q ) );
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

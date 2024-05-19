#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg


######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T}
     )::Int64  where {D,T}
    ratio::Float64 = 5.0;

    println( "|P| :", cardin( P ) );
    println( "|Q| :", cardin( Q ) );

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
            println( "Calling frechet_c_compute" );
            m = frechet_c_compute( P, Q );
            println( "... done" );
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
    println( "Reading files..." );
    P = Polygon_read_file( fl_a );
    Q = Polygon_read_file( fl_b );
    println( "Files read!..." );

    sgn = frechet_comp( P, Q,  );
end



num_args = length( ARGS );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

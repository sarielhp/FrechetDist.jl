#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables


######################################################################


function  test_files( base_dir, queries_file )
    df = CSV.read( queries_file, DataFrame, types=String, header=false );
    for  i in  1:nrow(df)
        s_a = df[i,1];
        s_b = df[i,2];

        fl_a = base_dir * df[i,1];
        fl_b = base_dir * df[i,2];
        rad = parse( Float64, df[i,3] );
        println( fl_a, " ", fl_b, " ", rad );

        poly_a = Polygon_read_file( fl_a );
        poly_b = Polygon_read_file( fl_b );

        m = frechet_c_approx( poly_a, poly_b, 5.0 );
        if  m.leash < rad
            println( "d(", s_a, ", ", s_b, ") < ", rad );
            continue;
        end
        lb = m.leash / m.ratio;
        if  lb >  rad
            println( "d(", s_a, ", ", s_b, ") > ", rad );
            continue;
        end
        println( "UNDECIDED" );

    end
    #print( df );
end



num_args = length( ARGS );

#p = Polygon_read_plt_file( "data/birds/1787_1.plt" );
#Polygon_write_to_file( p, "xx" );
#q = Polygon_read_plt_file(  "xx" );
#Polygon_write_to_file( q, "uu" );

#exit( -1 );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( -1 );
end

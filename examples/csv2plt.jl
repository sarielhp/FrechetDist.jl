#! julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon



######################################################################

function  plt_convert( filename, target )
    if  ! isfile( filename )
        println( "Error: No file ", filename );
        exit( -1 );
    end
    df=CSV.read( filename, DataFrame )

    poly = Polygon2F();
    for  i in 1:nrow(df)
        #println( typeof( df[i,1] ) );
        #println( typeof( df[i,2] ) );

        pnt = npoint( df[i,1],  df[i,2] );
        push!( poly, pnt );
    end
    write_to_file( poly, target );
    #println( p );
    #exit( -1 );
end


function (@main)(ARGS)
    num_args = length( ARGS );

    if   num_args == 0
        println( "csv2plt.jl [file1] [file2] ... " );
        exit( -1 );
    end

    for  i in 1:num_args
        target = replace( ARGS[ i ], ".csv" => ".plt" )

        println( "--- ", ARGS[ i ] );
        plt_convert( ARGS[ i ], target );
    end

    return  0;
end

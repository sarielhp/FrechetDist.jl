#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables


######################################################################

function  plt_convert( filename )
    if  ! isfile( filename )
        println( "Error: No file ", filename );
        exit( -1 );
    end
    s = readline( filename );
    if   s != "Geolife trajectory"
        return
    end
    p = Polygon_read_plt_orig_file( filename );
    println( "Rewriting: ", filename );
    Polygon_write_to_file( p, filename );
end



num_args = length( ARGS );

if   num_args == 0
    println( "plt2plt.jl [file1] [file2] ... " );
    exit( -1 );
end

for  i in 1:num_args
    plt_convert( ARGS[ i ] );
end

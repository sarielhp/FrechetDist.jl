#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using BasicCompGeometry
import BasicCompGeometry: write_plt
using CSV, DataFrames
using Downloads
using Dates, Printf;

function  convert( filename, out )
    df = CSV.read( filename, DataFrame );
    
    P = PntSeq2F();
    #@assert( length(x) == length(y) )
    for  i ∈ 1:nrow( df )
        if  ismissing( df.latitude[ i ] )
            continue;
        end
        if  ismissing( df.longitude[ i ] )
            continue;
        end
        x_p = Float64( df.latitude[ i ] );
        y_p = Float64( df.longitude[ i ] );
        push!( P, point( y_p, x_p ) );
    end
    #exit( -1 );
    write_plt( out, P );
    println( "Out length: ", length( P ) );
end

function (@main)(ARGS)
    count = 0;
    for i ∈ 2:length(ARGS)
        count = count + 1;
        outname = ARGS[1]*@sprintf( "out_%03d.plt", count );
       
        println( ARGS[ i ], "  ==>  ", outname )
        convert( ARGS[ i ], outname );
    end
    println( "bogi" );
end

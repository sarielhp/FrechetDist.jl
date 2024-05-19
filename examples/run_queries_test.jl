#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables


######################################################################

function frechet_decider( P::Polygon{D,T}, Q::Polygon{D,T},
                          r )::Int64  where {D,T}
    f_debug::Bool = true;

    if  Dist( first( P ), first( Q ) ) > r
        return  1;
    end
    if  Dist( last( P ), last( Q ) ) > r
        return  1;
    end
    ratio::Float64 = 5.0;

    for  i in 1:10
        f_debug &&  ( i > 5 )  &&   println( "ratio: ", ratio );
        m = frechet_c_approx( P, Q, ratio );
        if  m.leash < r
            return  -1;
        end
        lb = m.leash / m.ratio;
        if  lb > r
            return  1;
        end

        #    ratio = (r / lb min( m.ratio, 1.01 );
        ratio = ((r / lb) - 1.0) / 2.0 + 1.0; # min( m.ratio, 1.01 );
        ratio = min( ratio, 1.1 );
        if  ( ratio <= 1.001 )
            m = frechet_c_compute( P, Q );
            if  m.leash > r
                return  1;
            end
            if  m.leash < r
                return  -1;
            end
            return  0;
        end
#        if  i > 2
#            println( "RATIO  ", i, " : ", ratio );
#        end
    end
    println( "UNDECIDED" );

    return  0;
end


function  test_files( base_dir, queries_file )
    println( "QF: ", queries_file );
    df = CSV.read( queries_file, DataFrame, types=String, header=false );

    #    println( df );

    println( "Reading..." );
    PA = Vector{Polygon2F}();
    QA = Vector{Polygon2F}();
    rads = Vector{Float64}();

    nr = nrow(df);
#    nr = min( 10, nr );
    for  i in  1:nr
        s_a = df[i,1];
        s_b = df[i,2];

        fl_a = base_dir * df[i,1];
        fl_b = base_dir * df[i,2];
        rad = parse( Float64, df[i,3] );
        #println( fl_a, " ", fl_b, " ", rad );

        poly_a = Polygon_read_file( fl_a );
        poly_b = Polygon_read_file( fl_b );

        push!( PA, poly_a );
        push!( QA, poly_b );
        push!( rads, rad );
    end

    println( "Figuring out distances..." );
    for  i in  1:length( PA )
        println( base_dir * df[i,1], " vs ", base_dir * df[i,2] );
        flush( stdout );
        sgn = frechet_decider( PA[ i ], QA[ i ], rads[ i ] );
    end

    println( "Text completed on : ", queries_file );

    #print( df );
end



num_args = length( ARGS );


#exit( -1 );

function  test_files_from_file( filename, base_dir )
    println( "filename::: ", filename );
    lines = readlines( filename );
    nr = length( lines );
    #nr = min( 10, nr );
    for  i in 1:nr
        r = lines[ i ]
        test_files( base_dir, r );
    end

    println( "TEST COMPLETED SUCCESSFULLY!" );
end


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( -1 );
end

if   num_args == 3  &&  ( ARGS[ 1 ] == "file" )
    println( "ARGS[3]: ", ARGS[3 ] );
    test_files_from_file( ARGS[3], ARGS[2] );
    exit( -1 );
end
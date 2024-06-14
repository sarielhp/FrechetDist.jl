#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables

AtomicInt = Threads.Atomic{Int}

######################################################################

function frechet_decider( P::Polygon{D,T}, Q::Polygon{D,T},
                          r )::Int64  where {D,T}
    f_debug::Bool = false;

    f_debug  &&  println( "|P|: ", cardin( P ), " |Q|: ", cardin( Q ) );

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
        if  ( ratio <= 1.01 )
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
    f_debug  &&  println( "UNDECIDED" );
    @assert( false );

    return  0;
end


struct  PolygonsInDir
    polys::Vector{Polygon2F};
    d::Dict{String, Integer};
end;

function  Base.getindex( P::PolygonsInDir, s::String)
    return   P.polys[ P.d[ s ] ];
end

function   read_polygons_in_dir( base_dir )
    P = PolygonsInDir( Vector{Polygon2F}(), Dict( "empty" => -1 ) );
    for (root, dirs, files) in walkdir( base_dir )
        for file in files
            ( file == "dataset.txt" )  &&  continue;
            ( file == "mixoutALL_shifted.mat" )  && continue;
            ( file == "trajectories.names" )  && continue;
            
            println( "Reading: ", base_dir * file, "      \r" );
            poly = Polygon_read_file( base_dir * file );
            push!( P.polys, poly );
            P.d[ file ] = length( P.polys );
        end
    end
    return P;
end


function  test_files( PID, base_dir, queries_file, prefix,
                      count::AtomicInt,
                      i_second::Int64 = 1
                      )
    println( "test_file: ", queries_file, prefix );
    df = CSV.read( queries_file, DataFrame, types=String, header=false );

    #println( df );

    #println( "Reading..." );
    PA = Vector{Polygon2F}();
    QA = Vector{Polygon2F}();
    rads = Vector{Float64}();

    nr = nrow(df);

    #    nr = min( 10, nr );
    for  i in  1:nr
        #println( "i: ", i );
        s_a = df[i,1];
        s_b = df[i,2];

        fl_a = base_dir * df[i,1];
        fl_b = base_dir * df[i,2];
        #println( "=======================================" );
        rad = parse( Float64, df[i,3] );
        #println( fl_a, " ", fl_b, " ", rad );

        #poly_a = Polygon_read_file( fl_a );
        #poly_b = Polygon_read_file( fl_b );
        #println( "s_a : ", s_a );
        #println( "s_b : ", s_b );
        poly_a = PID[ s_a ];
        poly_b = PID[ s_b ];
        #println( "after?" );
        
        push!( PA, poly_a );
        push!( QA, poly_b );
        push!( rads, rad );
    end

    #println( "Figuring out distances..." );
    for  i in  1:length( PA )
        if  ( i < i_second )
            continue;
        end
        if  ( ( count[]  & 0xfff ) == 0xfff )
            println( count[], " T",
                Threads.threadid(),"D : ", prefix, i, "/", length( PA ), " ",
                base_dir * df[i,1], "   ", base_dir * df[i,2] );
            flush( stdout );
        end
        sgn = frechet_decider( PA[ i ], QA[ i ], rads[ i ] );

        Threads.atomic_add!( count, 1 )
    end

    println( "Text completed on : ", queries_file );

    #print( df );
end



num_args = length( ARGS );


#exit( -1 );

function  do_chunk( PID, lines, base_dir, nr,
                    count::AtomicInt,
                    i_second::Int64 = 1 )
#        test_files( base_dir, r, prefix, i_second );

    #nr = min( 10, nr );
    #parentindices( chunk )
    #println( "\n\n\nDO_CHUNK" );
    #println( lines );
    for  i in eachindex( lines )
        #println( "Before..." );
        i_orig = lines.offset1 + lines.stride1*i
        r = lines[ i ]
        prefix = @sprintf( "[%d/%d] ", i_orig, nr );
        test_files( PID, base_dir, r, prefix, count, i_second );
    end

    return  0;
end


function  test_files_from_file( filename, base_dir,
                                i_main::Int64 = 1,
                                i_second::Int64 = 1 )
    rlines = readlines( filename );

    println( "filename::: ", filename );
    println( "base_dir ", base_dir );

    PID = read_polygons_in_dir( base_dir );

    lines = rlines[i_main:length( rlines )];
    nr = length( lines );
    count = AtomicInt( 0 );

    if   ( i_second > 1 )
        prefix = @sprintf( "[%d/%d] ", i_main, nr );
        test_files( base_dir, lines[ i_main ], prefix, count, i_second );
        return;
    end


    nt = Threads.nthreads();
    #nt = 1;
    chunks = Iterators.partition(lines, length(lines) ÷ nt )
    tasks = map(chunks) do chunk
        #println( parentindices( chunk ) );
        Threads.@spawn do_chunk( PID, chunk, base_dir, nr, count, i_second );
    end;
    fetch.(tasks);

    #nr = min( 10, nr );
        #=
    for  i in 1:nr
        if  ( i > i_main )
            i_second = 1;
        end
        r = lines[ i ]
        prefix = @sprintf( "[%d/%d] ", i, nr );
        test_files( base_dir, r, prefix, i_second );
    end
        =#
    println( "TEST COMPLETED SUCCESSFULLY!" );
    flush( stdout );
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

if   num_args == 5  &&  ( ARGS[ 1 ] == "file" )
    println( "ARGS[3]: ", ARGS[3 ] );
    test_files_from_file( ARGS[3], ARGS[2], parse( Int64, ARGS[4] ),
                          parse( Int64, ARGS[ 5 ] ) );
    exit( -1 );
end

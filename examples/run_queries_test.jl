#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables
#using Profile
#using InteractiveUtils
#using ProfileView

AtomicInt = Threads.Atomic{Int}

######################################################################
###################################################################3

struct  PolygonsInDir
    PHA::Vector{PolygonHierarchy};
    polys::Vector{Polygon2F};
    widths::Vector{Float64};
    d::Dict{String, Integer};
    lock::ReentrantLock
    ph_lock::ReentrantLock
end;


function  GetIndex( P::PolygonsInDir, s::String)
    return   P.d[ s ];
end


function  Base.getindex( PID::PolygonsInDir, s::String)
    return   P.polys[ PID.d[ s ] ];
end

function  PID_init()
    return  PolygonsInDir( Vector{PolygonHierarchy}(),
                       Vector{Polygon2F}(), Vector{Float64}(),
        Dict( "empty" => -1 ), ReentrantLock(),
     ReentrantLock() );
end


function   PID_read_polygon( PID::PolygonsInDir,
                             filename::String,
                             key::String,
                             f_compute_hier::Bool = false )::Bool
    if  ! isfile( filename )
        println( "ERROR. Unable to open [", filename, "]" );
        return  false;
    end

    poly = Polygon_read_file( filename );
    width = frechet_width_approx( poly );

    if  ( cardin( poly ) <= 1 )
        println( "Bad poly file: \n\n", base_dir*file,"\n\n" );
    end
    @assert( cardin( poly ) > 1 );

    lock( PID.lock );
    push!( PID.polys, poly );
    push!( PID.widths, width );
    PID.d[ key ] = length( PID.polys );
    unlock( PID.lock );

    if  f_compute_hier
        ph = ph_init( poly, PID.ph_lock );
        lock( PID.lock );
        push!( PID.PHA,  ph );
        unlock( PID.lock );
    end

    return  true;
end


function  PID_init_simp_hierarchy( PID::PolygonsInDir, f_parallel::Bool )
    println( "\nInitialize simplification hierarchies..." );

    phx = ph_init( PID.polys[1], PID.ph_lock )
    for  i  in 1:length(PID.polys)
        lock( PID.lock )
        push!( PID.PHA, phx );
        unlock( PID.lock );
    end

    cnt_ph = AtomicInt( 0 );

    if  f_parallel
        println( "Doing this in parallel..." );
        Threads.@threads for  i  in 1:length(PID.polys)
            #println( i );
            poly = PID.polys[ i ];
            ph = ph_init( poly, PID.ph_lock )
            ph_compute_hierarchy( ph, 20.4,
                max( 100, round( Int64, cardin( poly ) / 4) ) );
            lock( PID.lock );
            PID.PHA[ i ] = ph;
            unlock( PID.lock );
            Threads.atomic_add!( cnt_ph, 1 )
            if  ( cnt_ph[] & 0xff) == 0
                print( cnt_ph[], "        \r" );
            end
        end
        return
    end

    for  i  in 1:length(PID.polys)
        poly = PID.polys[ i ];
        ph = ph_init( poly, PID.ph_lock )

        ph_compute_hierarchy( ph, 20.4,
            max( 100, round( Int64, cardin( poly ) / 4) ) );

        lock( PID.lock );
        PID.PHA[ i ] = ph;
        unlock( PID.lock );

        Threads.atomic_add!( cnt_ph, 1 )
        if  ( cnt_ph[] & 0xff) == 0
            print( cnt_ph[], "        \r" );
        end
    end
end

function   read_polygons_in_dir( base_dir, f_parallel::Bool )
    limit::Int64 = 200000;

    println( "Scanning for polygon files..." );
    list = Vector{String}();
    PID = PID_init();
    for (root, dirs, files) in walkdir( base_dir )
        for file in files
            ( file == "dataset.txt" )  &&  continue;
            ( file == "mixoutALL_shifted.mat" )  && continue;
            ( file == "trajectories.names" )  && continue;

            push!( list, file );
        end
    end

    count = AtomicInt( 0 );

    println( "Reading polygon files..." );
    if  ( f_parallel )
        Threads.@threads for  file  in list
            Threads.atomic_add!( count, 1 )
            if  ( count[] > limit  )
                break;
            end
            if  (count[] & 0xfff) == 0xfff
                println( count[], "  Reading: ", base_dir * file, "      \r" );
            end

            PID_read_polygon( PID, base_dir * file, file )
        end
    else
        for  file  in list
            Threads.atomic_add!( count, 1 )
            if  ( count[] > limit  )
                break;
            end
            if  (count[] & 0xfff) == 0xfff
                println( count[], "  Reading: ", base_dir * file, "      \r" );
            end

            PID_read_polygon( PID, base_dir * file, file )
        end
    end

    PID_init_simp_hierarchy( PID, f_parallel );

    return PID;
end


function  eq( a::Float64, b::Float64, tolerance::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (tolerance* (abs( a)  + abs(b) ))
end


##########################################################################33


function frechet_decider_PID_slow( PID, i, j, r )::Int64
    f_debug::Bool = false;

    P = PID.polys[ i ];
    Q = PID.polys[ j ];

    l_a =  Dist( first( P ), first( Q ) );
    if  l_a > r
        return  1;
    end
    l_b = Dist( last( P ), last( Q ) );
    if  l_b > r
        return  1;
    end

    ratio = 4.0;
    for  i in 1:100
        f_debug &&  ( i > 5 )  &&   println( "ratio: ", ratio );
        m = frechet_c_approx( P, Q, ratio );
        if  m.leash < r
            return  -1;
        end
        lb = m.leash / m.ratio;
        if  lb > r
            #println( "???????????" );
            return  +1;
        end

        ratio = ((r / lb) - 1.0) / 2.0 + 1.0; # min( m.ratio, 1.01 );
        ratio = min( ratio, 1.1 );
        if  ( ratio <= 1.01 )
            m = frechet_c_compute( P, Q );
            #println( "*** m.leash - r: ", m.leash - r );
            if  m.leash > r
                return  +1;
            end
            if  m.leash < r
                return  -1;
            end
            return  0;
        end
    end
    @assert( false );

    return  0;
end

function frechet_decider_PID( PID::PolygonsInDir, i::Int64,
                              j::Int64, r::Float64 )::Int64
    f_debug::Bool = false;
    f_verify::Bool = false;

    P = PID.polys[ i ];
    Q = PID.polys[ j ];

    #f_debug  &&  println( "\n\n@@@@@@@@a@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
    #f_debug  &&  println( "|P|: ", cardin( P ), " |Q|: ", cardin( Q ) );

    l_a =  Dist( first( P ), first( Q ) );
    ( l_a > r )  &&   return  1;
    lb = max( l_a, Dist( last( P ), last( Q ) ) );

    ( lb > r )   &&   return  1;

    P_ph = PID.PHA[ i ];
    Q_ph = PID.PHA[ j ];

    w_P = P_ph.widths[ 1 ];
    w_Q = Q_ph.widths[ 1 ];

    f_debug && print( "w_P  : ", w_P );
    f_debug && print( "    w_Q  : ", w_Q );

    ub = lb + w_P + w_Q;
    if  ub < r
        return  -1;
    end

    delta = min( abs( r - lb ), abs( r - ub ), (w_P + w_Q)/4.0 );# / 0.9;
    f_debug  &&  println( "Lower bound: ", lb, "\nUpper bound: ", ub,
                          "\nr: ", r );

    f_monotone::Bool = false;
    for  iters::Int64 in 1:14
        w_trg = delta / 2.0 #2.0 #1.5 # / 2.0
        PA, wP = ph_approx( P_ph, w_trg );
        QA, wQ = ph_approx( Q_ph, w_trg );

        if  f_monotone
            f_debug  &&  println( "MONOTONE" );
            m, PA_A, QA_A = frechet_mono_via_refinement_delta( PA, QA,
                                                               delta,
                                                            false );
            l_min = m.lower_bound;#leash / m.ratio;
            l_max = m.leash;
        else
            l_min, l_max = FEVER_compute_range( PA, QA, ub )
        end

        # the strange thing is that the max does not need to be equal...
        if  f_verify
            println( "VERIFY" );
            l_min_a, l_max_a = frechet_ve_r_compute_range( PA, QA,
                                                    20000000.0 + ub )
            if  ( ! eq( l_min, l_min_a, 0.00001 ) ) #  ||  ( l_max != l_max_a ) )
                println( "l_min  : ", l_min );
                println( "l_min_a: ", l_min_a );
                println( "l_max  : ", l_max );
                println( "l_max_a: ", l_max_a );
                Polygon_write_to_file( PA, "bad_a.txt" );
                Polygon_write_to_file( QA, "bad_b.txt" );
                @assert( l_min == l_min_a );
            end
        end

        f_debug && println( l_min, "...", l_max, "   r:", r );
        if  ! f_monotone
            if  ( ( iters > 0 )  &&  ( l_min < r < l_max )
                  &&  ( ( l_max - l_min ) > 4.0*delta ) )
                #delta = min( abs( l_min -r ), abs( l_max - r ), delta );
                f_debug  &&  println( "zoom in..." );
                f_monotone = true;
                continue;
            end
        end
        lb = max( lb, l_min - wP - wQ )
        ub = min( ub, l_max + wP + wQ )

        ( lb > r )  &&  return  +1;

        ( ub < r )  &&  return  -1;

        delta = min( abs( l_max - r ), abs( l_min - r ),
                     ub - r, r - lb,
                     delta / 3.0 );
    end

    println( "SHOGI!" );

    m = frechet_c_compute( P, Q )
    return  round(Int64, sign( m.leash - r ) )
end


function frechet_decider_PID_new( PID, i, j, r )::Int64
    f_debug::Bool = false;

    P = PID.polys[ i ];
    Q = PID.polys[ j ];


    #f_debug  &&
    #println( "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
    f_debug  &&  println( "|P|: ", cardin( P ), " |Q|: ", cardin( Q ) );

    l_a =  Dist( first( P ), first( Q ) );
    if  l_a > r
        return  1;
    end
    l_b = Dist( last( P ), last( Q ) );
    if  l_b > r
        return  1;
    end

    P_ph = PID.PHA[ i ];
    Q_ph = PID.PHA[ j ];

    local wP = P_ph.widths[ 1 ];
    local wQ = Q_ph.widths[ 1 ];

    f_debug  &&  print( "wP  : ", wP,  "  ", "    wQ  : ", wQ );

    dist = max( l_a, l_b );
    ub_start = ub = dist + wP + wQ;

    if  ( ub < r )
        return  -1;
    end

    lb = dist - wP - wQ;
    if  ( lb > r )
        return  +1;
    end

    ratio::Float64 = 5.0;
    delta = min( abs( r - lb ), abs( r - ub ), wP, wQ ) / 10.0;
    for  iters::Int64 in 1:200
        @assert( iters < 20 );
        w_trg = delta / 3.0 #1.5 # / 2.0
        PA, wP = ph_approx( P_ph, w_trg );
        QA, wQ = ph_approx( Q_ph, w_trg );

        if  ( wP > w_trg ) ||   ( wQ > w_trg )
            println( "wP: ", wP, " > ", w_trg );
            println( "wQ: ", wQ, " > ", w_trg );
        end
        @assert( wP <= w_trg );
        @assert( wQ <= w_trg );

        #println( "=========== ve_r " );
        @time l_min, l_max = frechet_ve_r_compute_range( PA, QA, 2000000.0+ub_start );
        println( "===========aprox" );
        @time  frechet_c_approx( PA, QA, 1.01 );


        if  ( l_min < r < l_max )
            delta = min( abs( l_min -r ), abs( l_max - r ), delta );
            println( "delta: ", delta );
            m, PA_A, QA_A = frechet_mono_via_refinement_delta( PA, QA,
                                                               delta / 3.0,
                                                            false );
            l_min = m.lower_bound;#leash / m.ratio;
            l_max = m.leash;
        end

        #println( iters, ": ", l_min, "...", l_max, "  r: ", r );

        if  f_debug
            a_dist = frechet_c_compute( P, PA );
            b_dist = frechet_c_compute( PA, QA );
            c_dist = frechet_c_compute( QA, Q );
            all_dist = frechet_c_compute( P, Q, false );
            println( "dist( P, PA ): ", a_dist.leash );
            println( "dist( PA, QA ): ", b_dist.leash );
            println( "dist( QA, Q ): ", c_dist.leash );
            println( "dist( P, Q ): ", all_dist.leash );
            m_a = Morphing_combine( a_dist, b_dist );
            m_b = Morphing_combine( m_a, c_dist );
            println( "merged leash: ", m_b.leash );
            m_new = frechet_mono_via_refinement( P, Q, 1.0001 )[1]
            println( "d(P,Q) via refinement  : ", m_new.leash );

            m_verify = frechet_ve_r_compute( PA, QA );
            l_min_2 = m_verify.leash;
            m_m_verify = Morphing_monotonize( m_verify )
            l_max_2 = m_m_verify.leash;
            println( l_min_2, "...", l_max_2, "   L: ", l_max_2 - l_min_2 );
        end

        lb = l_min - wP - wQ
        ub = l_max + wP + wQ

        #println( lb, "...", ub );
        if  ( lb > r )
            return  +1;
        end
        if  ( ub < r )
            return  -1;
        end

        f_debug  &&  println( "ub B : ", ub );

        delta = min( min( abs( l_min - r ), abs( l_max - r ),
                          abs( lb - r ),    abs( ub - r ) ), delta )/ 4.0 ;
    end

    #=
    println( "SHOGI!" );

    for  i in 1:10
        f_debug &&  ( i > 5 )  &&   println( "ratio: ", ratio );
        m = frechet_c_approx( P, Q, ratio );

        ( m.leash < r )  &&  return  -1;

        lb = m.leash / m.ratio;

        ( lb > r )  &&  return  +1;

        ratio = ((r / lb) - 1.0) / 2.0 + 1.0;
        ratio = min( ratio, 1.1 );
        if  ( ratio <= 1.01 )
            m = frechet_c_compute( P, Q );
            return  round( Int64, sign( m.leash - r ) );
        end
    end
    =#
    f_debug  &&  println( "UNDECIDED" );
    @assert( false );

    return  0;
end


function frechet_decider( P::Polygon{D,T}, Q::Polygon{D,T},
    r, w_P::Float64 = -1.0, w_Q::Float64 = -1.0,
    f_use_widths::Bool = false;
)::Int64  where {D,T}
    f_debug::Bool = false;

    f_debug  &&  println( "|P|: ", cardin( P ), " |Q|: ", cardin( Q ) );

    l_a =  Dist( first( P ), first( Q ) );
    if  l_a > r
        return  1;
    end
    l_b = Dist( last( P ), last( Q ) );
    if  l_b > r
        return  1;
    end
    if  f_use_widths
        ub = max( l_a, l_b ) + w_P + w_Q;
        if  ub < r
            return  -1;
        end
        lb = max( l_a, l_b ) - w_P - w_Q;
        if  lb > r
            return  +1;
        end
    end

    ratio::Float64 = 5.0;

    for  i in 1:10
        println( "Iter: ", i );
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
    end
    f_debug  &&  println( "UNDECIDED" );
    @assert( false );

    return  0;
end

mutable struct  test_info_t
    i_P::Int64
    i_Q::Int64
    rad::Float64;
    f_l_a::String
    f_l_b::String
    runtime::Float64
end

function Base.isless(a::test_info_t, b::test_info_t)
    return  isless(b.runtime, a.runtime)
end


function  run_tests( PID::PolygonsInDir, tests::Vector{test_info_t},
    count::AtomicInt, f_verify::Bool,
    rng::UnitRange{Int64} = 0:0,
)
    errors::Int64 = 0;

    mask::Int64 = 0x1fff; #0x1; #0x1fff
    #println( "run_tests..." );

    if  ( rng == 0:0 )
        rng = 1:length( tests );
    end

    f_verify  &&  println( "Verification is on, test would run much slower..." );
    for  i in  rng
        #=if  i > 10000000
            return 0;
        end=#
        if  ( ( count[]  & mask ) == mask )
            println( count[], " T",
                Threads.threadid(),"D : ", i, "/", length( tests ) );
            flush( stdout );
        end
        t = tests[ i ];
        #println( "Before frechet_decider..." );
        #println( i,":   fl_a: ", t.f_l_a, "   fl_b: ", t.f_l_b );
        ms = @timed sgn = frechet_decider_PID( PID, t.i_P, t.i_Q, t.rad );
        t.runtime = ms.time;

        if  ( f_verify )
            #println( "verifying!" );
            sgn_slow = frechet_decider_PID_slow( PID, t.i_P, t.i_Q, t.rad )
            if  ( sgn != sgn_slow )
                P = PID.polys[ t.i_P ];
                Q = PID.polys[ t.i_Q ];
                m = frechet_c_compute( P, Q );

                sgn_real::Int64 = round(Int64, sign( m.leash - t.rad ) );
                if  ( sgn != sgn_real )
                    errors = errors + 1;
                    println( "sgn      : ", sgn );
                    println( "sgn_slow : ", sgn_slow );
                    println( "sgn_real : ", sgn_real );
                    println( "r        : ", t.rad );
                    println( "m.leash  : ", m.leash );
                    #@assert( sgn == sgn_slow );
                    open( "errors_log.txt", "a" ) do file
                        println( file, t.f_l_a, " ", t.f_l_b, " ", t.rad );
                    end
                end
            end
        end

        Threads.atomic_add!( count, 1 )
        #if ( count[] > 1000 )
        #    return;
        #end;
    end

    sort!( tests );

    open( "tests.txt", "w" ) do file
        for  t in tests
            println( file, t.f_l_a, " ", t.f_l_b, " ", t.rad );
        end
    end
    open( "tests_times.txt", "w" ) do file
        for  t in tests
            println( file, t.f_l_a, " ", t.f_l_b, " ", t.rad, " ", t.runtime );
        end
    end

    return  errors;
end

function  read_query_file( base_dir, queries_file, tests, PID )
    df = CSV.read( queries_file, DataFrame, types=String, header=false );

    nr = nrow(df);

    for  i in  1:nr
        s_a = df[i,1];
        s_b = df[i,2];

        fl_a = base_dir * df[i,1];
        fl_b = base_dir * df[i,2];

        rad = parse( Float64, df[i,3] );
        if  ! haskey( PID.d, s_a )
            if  ( !PID_read_polygon( PID, fl_a, s_a, true ) )
                continue;
            end
        end
        if  ! haskey( PID.d, s_b )
            if  ( !PID_read_polygon( PID, fl_b, s_b, true ) )
                continue;
            end
        end
        ind_a = GetIndex( PID, s_a );
        ind_b = GetIndex( PID, s_b );

        t = test_info_t( ind_a, ind_b, rad, fl_a, fl_b, 0.0 );
        push!( tests, t );
    end
end


function  test_files( PID, base_dir, queries_file, prefix,
                      count::AtomicInt,
                      f_verify::Bool )::Int64
    errors::Int64 = 0;

    println( prefix, " : ", queries_file );

    tests = Vector{test_info_t}();
    #println( "Reading query file..." );
    read_query_file( base_dir, queries_file, tests, PID )

    errors = run_tests( PID, tests, count, f_verify )
    #@profile    errors = run_tests( PID, tests, count, f_verify )
    #Profile.print()

    if  ( errors > 0 )
        println( "ERRORS: ", errors );
    end

    return  errors;
end

function  test_single_file( filename, f_verify::Bool )
    tests = Vector{test_info_t}();
    PID = PID_init();

    println( "Reading query file..." );
    read_query_file( "", filename, tests, PID )

    count = AtomicInt( 0 );

    return  run_tests( PID, tests, count, false );
end





function  do_array( PID, lines, base_dir, nr,
                    count::AtomicInt,
                    f_verify::Bool )
    errors::Int64 = 0;

    tests = Vector{test_info_t}();

    println( "\nReading query files..." );
    for  line  in  lines
        read_query_file( base_dir, line, tests, PID )
    end

    ### Do the tests...
    println( "\nRunning tests..." );
    errors = run_tests( PID, tests, count, f_verify )

    ( errors > 0 )  &&  println( "ERRORS TOTAL: ", errors );

end

function  do_chunk( PID, lines, base_dir, nr,
                    count::AtomicInt )
    for  i in eachindex( lines )
        #println( "Before..." );
        i_orig = lines.offset1 + lines.stride1*i
        r = lines[ i ]
        prefix = @sprintf( "[%d/%d] ", i_orig, nr );
        test_files( PID, base_dir, r, prefix, count, false );
    end

    return  0;
end


function  test_files_from_file( filename, base_dir,
    f_verify::Bool,
    f_serial::Bool = false
)
    rlines = readlines( filename );

    println( "\nReading polygons..." );
    PID = read_polygons_in_dir( base_dir, ! f_serial );

    println( "\n" );
    lines = rlines[1:length( rlines )];
    nr = length( lines );
    count = AtomicInt( 0 );

    if  ( f_serial )
        @time do_array( PID, lines, base_dir, nr, count, f_verify )
    else
        nt = Threads.nthreads();
        chunks = Iterators.partition(lines, length(lines) ÷ nt )
        tasks = map(chunks) do chunk
            #println( parentindices( chunk ) );
            Threads.@spawn do_chunk( PID, chunk, base_dir, nr, count );
        end;
        fetch.(tasks);
    end

    println( "TEST COMPLETED SUCCESSFULLY!" );
    println( "# of pairs compared : ", count[] );
    flush( stdout );
end

########################################################################
########### Main

f_verify_run::Bool = false;

num_args = length( ARGS );

if   num_args == 2  &&  ( ARGS[ 1 ] == "test_file" )
    test_single_file( ARGS[2], f_verify_run );
    exit( 0 );
end
if   num_args == 3  &&  ( ARGS[ 1 ] == "file" )
    test_files_from_file( ARGS[3], ARGS[2], f_verify_run );
    exit( 0 );
end

if   num_args == 3  &&  ( ARGS[ 1 ] == "sfile" )
    test_files_from_file( ARGS[3], ARGS[2], f_verify_run, true );
    exit( 0 );
end

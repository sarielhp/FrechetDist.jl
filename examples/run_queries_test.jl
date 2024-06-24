#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables
#using Profile

AtomicInt = Threads.Atomic{Int}

######################################################################


mutable struct  PolygonHierarchy
    len::Float64;
    P::Polygon2F;
    plt::Vector{Float64};

    polys::Vector{Polygon2F};
    widths::Vector{Float64};

end

function ph_push!( ph::PolygonHierarchy, Q::Polygon2F, w::Float64 )
    @assert( cardin( Q ) > 1 )
    push!( ph.polys, Q );
    push!( ph.widths, w );
end


function  ph_push_target_exp( ph::PolygonHierarchy, w::Float64,
                              P::Polygon2F, lmt::Int64, w_extra::Float64 )
    T, T_indices = frechet_simplify_w_exp( P, w );
    if  ( cardin( T ) > lmt )
        return  -1.0;
    end

    mT = frechet_c_mono_approx_subcurve( P, T, T_indices )[ 1 ];
    new_w = mT.leash + w_extra;

    ph_push!( ph, T, new_w );

    return  new_w;
end


function  ph_push_target( ph::PolygonHierarchy, w::Float64 )
    P = ph.P;

    if  ( cardin( P ) <= 10 )
        return  true;
    end

    R, R_indices = frechet_simplify_to_width( P, w );
    ###T, T_indices = frechet_simplify_w_exp( P, w );

    # No point continuing...
    if  ( ( 1.1 * cardin( R ) ) > cardin( P ) )
        return  true;
    end
    mr = frechet_c_mono_approx_subcurve( P, R, R_indices )[ 1 ];
    #println( "## :", cardin( P ), "  # :", cardin( R ), "    R_w: ", mr.leash );

    push!( ph.polys, R );
    push!( ph.widths, mr.leash );

    return  false;
end


function  ph_print( ph::PolygonHierarchy )
    P = ph.P
    println( "--------------------------------------------" );
    println( "|P|: ", cardin( P ) );
    for  i  in  1:length( ph.widths )
        println( "w[", i, "] :", ph.widths[ i ],
                 "  #: ", cardin( ph.polys[ i ] ) );
    end
    println( "----------------------------------------------\n" );
end


function   ph_approx( ph::PolygonHierarchy, w::Float64 )

    resolution::Float64 = 4.0;
    simp_threshold::Float64 = 4.0;

    #resolution, simp_rsolution = 1.1, 4.0; # 242 seconds
    # resolution, simp_rsolution = 4.0, 4.0; # 205
    #resolution, simp_rsolution = 8.0, 4.0; # 213
    #resolution, simp_rsolution = 1.4, 4.0; # 215
    #resolution, simp_rsolution = 1.4, 2.0; # 215
    resolution, simp_rsolution = 1.4, 10.0; # 214

    #println( "---------------------" );
    for  i  in  1:length(ph.widths)
        #println( ph.widths[ i ] );
        if  ( ph.widths[ i ] < w <= ( resolution * ph.widths[ i ] ) )
            #println( "Bingo!" );
            return  ph.polys[ i ], ph.widths[ i ]
        end
    end

    Z, Z_indices, wZ = frechet_approx_from_pallete( ph.P, ph.plt,
                                                    w / simp_threshold );
    #w = ph_push_target_exp( phA,   max( wA - wZ, 0.0 ), Z, lmt, wZ )

    wtmp = max( w-wZ, 0.0 );

    X, X_indices = frechet_simplify_w_exp( Z, wtmp );
    mX = frechet_c_mono_approx_subcurve( Z, X, X_indices )[ 1 ];

    w_out = mX.leash + wZ;

    push!( ph.polys, X );
    push!( ph.widths, w_out );

    return  X, mX.leash + wZ;
end


function ph_init( P::Polygon2F )
    ph = PolygonHierarchy( 0.0, P, Vector{Int64}(), Vector{Polygon2F}(),
                           Vector{Float64}() );
    ph.len = Polygon_length( P );
    ph.plt = frechet_pallete( P );

    w = frechet_width_approx( P );
    ph_push!( ph, Polygon_spine( P ), w );

    return  ph;
end


#global mx::Float64 = 0.0;

function  compute_simp_hierarchy( P::Polygon2F )
    #ph = ph_init( P );
    phA = ph_init( P );

    w = phA.widths[ 1 ];

    #wL = 0.0;
#    PL = P;
    #=
    ratio::Float64 = 1.4
    lmt::Int64 = min( round( Int64, cardin( P ) * 0.95 ), 500 );
    #=f_large::Bool = false;
    if  cardin( P ) > 1000
        #println( "LARGE" );
        ratio = 8.0;
        lmt = 500;
        wL = w / 1000.0;
        PL, PL_indices = frechet_simplify_to_width( P, wL );
        mr = frechet_c_mono_approx_subcurve( P, PL, PL_indices )[ 1 ];
        wL = mr.leash;
        #println( "|PL| ", cardin( PL ), "  |P|: ", cardin( P ),
        #         "  wL: ", wL );
        f_large = true;
    elseif  cardin( P ) > 500
        reatio = 2.0;
        lmt = 250;
    elseif  cardin( P ) > 100
        lmt = 50;
    end=#

    for  i in 1:200
        if cardin( last( phA.polys ) ) >= lmt
            break
        end
        wA = w / ratio;

        Z, Z_indices, wZ = frechet_approx_from_pallete( P, phA.plt, wA / 10.0 );

        w = ph_push_target_exp( phA,   max( wA - wZ, 0.0 ), Z, lmt, wZ )
#        else
#            w = ph_push_target_exp( phA,     wA, P, lmt, 0.0 )
#        end
        ( w < -0.5 )  &&  break;
    end
#    ph_push!( phA, P, 0.0 );
    #ph_print( phA );

    #f_large && ph_print( P, phA );
    =#

    return  phA;
end


struct  PolygonsInDir
    PHA::Vector{PolygonHierarchy};
    polys::Vector{Polygon2F};
    widths::Vector{Float64};
    d::Dict{String, Integer};
end;


function  GetIndex( P::PolygonsInDir, s::String)
    return   P.d[ s ];
end


function  Base.getindex( P::PolygonsInDir, s::String)
    return   P.polys[ P.d[ s ] ];
end

function   read_polygons_in_dir( base_dir, f_parallel::Bool )
    limit::Int64 = 500000;

    count::Int64 = 0;
    P = PolygonsInDir( Vector{PolygonHierarchy}(),
                       Vector{Polygon2F}(), Vector{Float64}(),
                       Dict( "empty" => -1 ) );
    for (root, dirs, files) in walkdir( base_dir )
        if  ( count > limit )
            break;
        end
        for file in files
            ( file == "dataset.txt" )  &&  continue;
            ( file == "mixoutALL_shifted.mat" )  && continue;
            ( file == "trajectories.names" )  && continue;

            count = count + 1;
            if  ( count > limit  )
                break;
            end
            if  (count & 0xfff) == 0xfff
                println( count, "  Reading: ", base_dir * file, "      \r" );
            end
            poly = Polygon_read_file( base_dir * file );
            width = frechet_width_approx( poly );

            @assert( cardin( poly ) > 1 );

            push!( P.polys, poly );
            push!( P.widths, width );
            P.d[ file ] = length( P.polys );
        end
    end

    phx = compute_simp_hierarchy( P.polys[1] )
    for  i  in 1:length(P.polys)
        push!( P.PHA, phx );
    end

    println( "Computing simplification hierarchies..." );
    cnt_ph = AtomicInt( 0 );

    if  f_parallel
        Threads.@threads for  i  in 1:length(P.polys)
            #println( i );
            poly = P.polys[ i ];
            ph = compute_simp_hierarchy( poly )
            P.PHA[ i ] = ph;
            Threads.atomic_add!( cnt_ph, 1 )
            if  ( cnt_ph[] & 0xff) == 0
                print( cnt_ph[], "        \r" );
            end
        end
    else
        for  i  in 1:length(P.polys)
            #println( i );
            poly = P.polys[ i ];
            ph = compute_simp_hierarchy( poly )
            P.PHA[ i ] = ph;
            Threads.atomic_add!( cnt_ph, 1 )
            if  ( cnt_ph[] & 0xff) == 0
                print( cnt_ph[], "        \r" );
            end
        end
    end
    return P;
end


function  eq( a::Float64, b::Float64, tolerance::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (tolerance* (abs( a)  + abs(b) ))
end


##########################################################################33


function frechet_decider_PID_old( PID, i, j, r )::Int64
    f_debug::Bool = false;

    P = PID.polys[ i ];
    Q = PID.polys[ j ];

    f_debug && println( "\n\n@@@@@@@@a@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
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

    w_P = P_ph.widths[ 1 ];
    w_Q = Q_ph.widths[ 1 ];

    f_debug && print( "w_P  : ", w_P );
    f_debug && print( "    w_Q  : ", w_Q );

    ub = max( l_a, l_b ) + w_P + w_Q;
    if  ub < r
        return  -1;
    end
    lb = max( l_a, l_b ) - w_P - w_Q;
    if  lb > r
        return  +1;
    end

    ratio::Float64 = 5.0;
    delta = min( abs( r - lb ), abs( r - ub ) );
    mi = max( length( P_ph.polys ), length( Q_ph.polys ) );
    for  iv in 2:mi
        iP = min( iv, length( P_ph.widths ) );
        iQ = min( iv, length( Q_ph.widths ) );
        w_P = P_ph.widths[ iP ];
        w_Q = Q_ph.widths[ iQ ];

        if  ( ( w_P + w_Q ) > delta )
            continue;
        end
        P_a = P_ph.polys[ iP ];
        Q_a = Q_ph.polys[ iQ ];

        m_leash = frechet_ve_r_compute_mono_dist( P_a, Q_a, ub );

        #=
        m = frechet_ve_r_compute( P_a, Q_a );
        mm_leash =  Morphing_monotone_leash( m );

        if  ( ! eq( m_leash, mm_leash, 0.000001 ) )
            println( "XXX :", m_leash, "   ", mm_leash );
            exit( -1 );
        end
        =#
        ## m_leash = frechet_ve_r_compute( P_a, Q_a );
        lb = m_leash - w_P - w_Q


        if  f_debug
            println( "---------------------------------------------------" );
            println( "|P_a|: ", cardin( P_a ) );
            println( "|Q_a|: ", cardin( Q_a ) );
            println( "r    : ", r );
            println( "ve_l : ", m_leash );
            println( "w_P  : ", w_P );
            println( "w_Q  : ", w_Q );
            println( "lb A : ", lb );
        end
        if  ( lb > r )
            return  +1;
        end
        lb = m_leash - w_P - w_Q
        f_debug  &&  println( "lb B : ", lb );
        if  ( lb > r )
            return  +1;
        end

        ub = m_leash + w_P + w_Q
        f_debug  &&  println( "ub B : ", ub );
        if  ( ub < r )
            return  -1;
        end
        delta = abs( m_leash - r );
    end

#    println( "SHOGI!" );

    for  i in 1:10
#        println( "Iter: ", i );
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

function frechet_decider_PID( PID, i, j, r )::Int64
    f_debug::Bool = false;

    P = PID.polys[ i ];
    Q = PID.polys[ j ];

    f_debug && println( "\n\n@@@@@@@@a@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
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

    w_P = P_ph.widths[ 1 ];
    w_Q = Q_ph.widths[ 1 ];

    f_debug && print( "w_P  : ", w_P );
    f_debug && print( "    w_Q  : ", w_Q );

    ub = max( l_a, l_b ) + w_P + w_Q;
    if  ub < r
        return  -1;
    end
    lb = max( l_a, l_b ) - w_P - w_Q;
    if  lb > r
        return  +1;
    end

    ratio::Float64 = 5.0;
    delta = min( abs( r - lb ), abs( r - ub ) );
    #mi = max( length( P_ph.polys ), length( Q_ph.polys ) );
    println( "\n" );
    while  true
        w_trg = delta / 2.0 #1.5 # / 2.0
        PA, wP = ph_approx( P_ph, w_trg );
        QA, wQ = ph_approx( Q_ph, w_trg );

        print( "delta: ", delta, "  wP: ", wP, "  wQ: ", wQ, "  " );
        println( cardin( PA ), " / ", cardin( P ),  "   |   ",
                 cardin( QA ), " / ", cardin( Q ) );

        m_leash = frechet_ve_r_compute_mono_dist( PA, QA, ub );

        lb = m_leash - wP - wQ


        if  f_debug
            println( "---------------------------------------------------" );
            println( "|P_a|: ", cardin( PA ) );
            println( "|Q_a|: ", cardin( QA ) );
            println( "r    : ", r );
            println( "ve_l : ", m_leash );
            println( "w_P  : ", wP );
            println( "w_Q  : ", wQ );
            println( "lb A : ", lb );
        end
        if  ( lb > r )
            return  +1;
        end

        ub = m_leash + wP + wQ
        f_debug  &&  println( "ub B : ", ub );
        if  ( ub < r )
            return  -1;
        end
        delta = abs( m_leash - r );
    end

#    println( "SHOGI!" );

    for  i in 1:10
#        println( "Iter: ", i );
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
#        if  i > 2
#            println( "RATIO  ", i, " : ", ratio );
#        end
    end
    f_debug  &&  println( "UNDECIDED" );
    @assert( false );

    return  0;
end

struct  test_info_t
    i_P::Int64
    i_Q::Int64
    rad::Float64;
end

function  test_files( PID, base_dir, queries_file, prefix,
                      count::AtomicInt,
                      i_second::Int64 = 1
                      )
    println( "test_file: ", queries_file, prefix );
    df = CSV.read( queries_file, DataFrame, types=String, header=false );

    tests = Vector{test_info_t}();

    nr = nrow(df);

    #    nr = min( 10, nr );
    for  i in  1:nr
        #println( "i: ", i );
        s_a = df[i,1];
        s_b = df[i,2];

        fl_a = base_dir * df[i,1];
        fl_b = base_dir * df[i,2];
        # println( fl_a, "  ", fl_b );
        #println( "=======================================" );
        rad = parse( Float64, df[i,3] );
        if  ! haskey( PID.d, s_a )
            continue;
        end
        if  ! haskey( PID.d, s_b )
            continue;
        end
        ind_a = GetIndex( PID, s_a );
        ind_b = GetIndex( PID, s_b );

        t = test_info_t( ind_a, ind_b, rad );
        push!( tests, t );
    end

    for  i in  1:length( tests )
        if  ( i < i_second )
            continue;
        end
        if  ( ( count[]  & 0x1fff ) == 0x1fff )
            println( count[], " T",
                Threads.threadid(),"D : ", prefix, i, "/", length( tests ), " ",
                base_dir * df[i,1], "   ", base_dir * df[i,2] );
            flush( stdout );
        end
        t = tests[ i ];
        sgn = frechet_decider_PID( PID, t.i_P, t.i_Q, t.rad );

        Threads.atomic_add!( count, 1 )
        #if ( count[] > 1000 )
        #    return;
        #end;
    end

    #println( "Text completed on : ", queries_file );

    #print( df );
end



num_args = length( ARGS );



function  do_array( PID, lines, base_dir, nr,
                    count::AtomicInt,
                    i_second::Int64 = 1 )
    for  i in eachindex( lines )
        r = lines[ i ]
        prefix = @sprintf( "[%d/%d] ", i, length( lines ) );
        test_files( PID, base_dir, r, prefix, count, i_second );
        #if ( count[] > 1000 )
        #    return;
        #end;
    end
end

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
    i_second::Int64 = 1,
    f_serial::Bool = false
)
    rlines = readlines( filename );

    println( "filename::: ", filename );
    println( "base_dir ", base_dir );

    PID = read_polygons_in_dir( base_dir, ! f_serial );

    lines = rlines[i_main:length( rlines )];
    nr = length( lines );
    count = AtomicInt( 0 );

    if   ( i_second > 1 )
        prefix = @sprintf( "[%d/%d] ", i_main, nr );
        test_files( base_dir, lines[ i_main ], prefix, count, i_second );
        return;
    end


    if  ( f_serial )
        @time do_array( PID, lines, base_dir, nr, count, i_second )
    else
        nt = Threads.nthreads();
        chunks = Iterators.partition(lines, length(lines) ÷ nt )
        tasks = map(chunks) do chunk
            #println( parentindices( chunk ) );
            Threads.@spawn do_chunk( PID, chunk, base_dir, nr, count, i_second );
        end;
        fetch.(tasks);
    end

    println( "TEST COMPLETED SUCCESSFULLY!" );
    println( "# of pairs compared : ", count[] );
    flush( stdout );
end


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( -1 );
end

if   num_args == 3  &&  ( ARGS[ 1 ] == "file" )
    println( "ARGS[3]: ", ARGS[3 ] );
    test_files_from_file( ARGS[3], ARGS[2] );
    exit( 0 );
end

if   num_args == 3  &&  ( ARGS[ 1 ] == "sfile" )
    println( "ARGS[3]: ", ARGS[3 ] );
    test_files_from_file( ARGS[3], ARGS[2], 1, 1, true );
    exit( 0 );
end

if   num_args == 5  &&  ( ARGS[ 1 ] == "file" )
    println( "ARGS[3]: ", ARGS[3 ] );
    test_files_from_file( ARGS[3], ARGS[2], parse( Int64, ARGS[4] ),
                          parse( Int64, ARGS[ 5 ] ) );
    exit( -1 );
end

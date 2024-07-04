##################################################################3
# fever.jl
#
# An attempt to implement a faster version of the VE Retractable
# Fréchet distance. The idea is to use static arrays to minimize
# dynamic algorithms.
##################################################################3

using DataStructures

const ID_START::Int64 = 1
const ID_END::Int64 = 2

struct EIDCalc
    nP::Int64
    nQ::Int64
end

function   ID_init( nP::Int64, nQ::Int64 )
    e = EIDCalc( nP, nQ );

    return  e;
end

function   ID_val_inner( i::Int64, j::Int64, r::Int64, nP::Int64 )
    @assert( ( r == 0 )  ||  (r == 1 ) )
    return (( j * nP + i ) << 1)  |  r;
end

function   ID_get( e::EIDCalc, i::Int64 = 0, i_is_vert::Bool = false,
    j::Int64 = 0, j_is_vert::Bool = false )
    if  ( i_is_vert  &&  j_is_vert )
        ( ( i == 1 )  &&  ( j == 1 ) )  &&  return ID_START;
        ( ( i == e.nP )  &&  ( j == e.nQ ) )  &&  return ID_END;
        @assert( false );
    end
    @assert( i_is_vert ||  j_is_vert );

    r::Int64 = ( i_is_vert ) ? 0 : 1 ;
    # i ∈ 1...nP-1, j ∈ 1.. nQ-1

    return  ID_val_inner( i, j, r, e.nP );
end

mutable  struct FeverCoords
    i::Int64;
    i_is_vert::Bool;
    j::Int64;
    j_is_vert::Bool;
end

function   ID_get_fields( e::EIDCalc, id::Int64, fc::FeverCoords )::Int64
    if  ( id == ID_START)
        fc.i = 1;
        fc.j = 1;
        fc.i_is_vert = true;
        fc.j_is_vert = true;
        return 0
    elseif  ( id == ID_END)
        fc.i = e.nP;
        fc.j = e.nQ;
        fc.i_is_vert = true;
        fc.j_is_vert = true;
        return 0
    end
    fc.i_is_vert::Bool = ( ( id & 1 ) == 0 )
    fc.j_is_vert = ! fc.i_is_vert;
    ix = id >> 1;
    fc.j = div( (ix - 1), e.nP );
    fc.i = ix - fc.j * e.nP;
    return  0;
end

function  ID_get_max( e::EIDCalc )
    #println( e.nP, e.nQ );
    return   ID_get( e, e.nP - 1, false, e.nQ - 1, true ) + 4;
end

function  ID_tester( nP::Int64, nQ::Int64)
    e = ID_init( nP, nQ );
    v = Vector{Int64}()

    mx = ID_get_max( e );
    push!( v, ID_get( e, 1, true, 1, true ) );
    push!( v, ID_get( e, nP, true, nQ, true ) );

    for  i in 1:(nP-1)
        for  j in 1:(nQ-1)
            id = ID_get( e, i, false, j, true );
            push!( v, ID_get( e, i, false, j, true ) );
            if  ( ID_status( id ) == 0 )
                println( "i: ", i ) ;
                println( "j: ", j ) ;
                println( "id: ", id );
                println( "id status: ", ID_status( id ) );

                @assert( ID_status( id ) == 1 );
            end
            fc = FeverCoords( 0, false, 0, false );
            ID_get_fields( e, last( v ), fc );
            @assert( i == fc.i );
            @assert( j == fc.j );
            @assert( fc.i_is_vert == false );
            @assert( fc.j_is_vert == true );
            
            push!( v, ID_get( e, i, true, j, false ) );
            if  ( ID_status( last( v  ) ) == 1 )
                println( "i: ", i ) ;
                println( "j: ", j ) ;
                println( "id: ", last( v ) );
                println( "id status: ", ID_status( last( v ) ) );

                @assert( ID_status( last(v) ) == 0 );
            end
            ID_get_fields( e, last( v ), fc );
            @assert( i == fc.i );
            @assert( j == fc.j );
            @assert( fc.i_is_vert == true );
            @assert( fc.j_is_vert == false );
        end
    end

    sort!( v );
    @assert( 1 <= v[ 1 ] <= mx );
    for  i  in 2:length( v )
        @assert( 1 <= v[ i ] <= mx );
        if  v[ i - 1 ] == v[ i ]
            println( "ERROR!" );
            @assert( v[ i - 1 ] != v[ i ] );
        end
    end
    #println( "nP: ", nP, " nQ: ", nQ, "LB: ", nP * nQ * 2, " vs. ", length( v ), "  mx: ", mx );

    # DEBUG
    #println( "Test passed" );
    return  0;
end

mutable struct  EventsOrder{T} <: Base.Order.Ordering
    values::Vector{T}
end

import Base.Order.lt
function  lt(o::EventsOrder, a::Int64, b::Int64 )
    return  isless( o.values[ a ], o.values[ b ] );
end

mutable struct FEVERContext{N,T}
    P::Polygon{N,T};
    Q::Polygon{N,T};

    eid::EIDCalc;
    vals::Vector{T}   # An array containing all the event values.
    events_queue::Vector{Int64};  # an array containing the IDs of all
                                  # the events in the heap
    comparator::EventsOrder{T};
    heap::BinaryHeap{Int64,EventsOrder{T}};
    #::BinaryHeap{EventsOrder{T},Int64};

    prev::Vector{Int64};

    n_p::Int64 #0
    n_q::Int64 #0
    f_upper_bound::Bool
    upper_bound::T
end

function   FEVER_Context(  P::Polygon{N,T},
    Q::Polygon{N,T} ) where {N,T}

    @assert( ( cardin( P ) > 1 )  &&  ( cardin( Q ) > 1 ) );

    #println( "FEVER_context start..." );
    eid = ID_init( cardin( P ), cardin( Q ) );
    mx = ID_get_max( eid );
    #@time println( "Zeros..." );
    #vals = zeros( T, mx );
    vals = fill( T( -1.0 ), mx );
    prev = zeros( Int64, mx );


    vals[ ID_START ] = Dist( first( P ), first( Q ) );
    vals[ ID_END ] = Dist( last( P ), last( Q ) );

    #println( "\nWWW" );
    comparator = EventsOrder{T}( vals );

    events_queue = Vector{Int64}();
    #println( "Heap init..." );
    heap = BinaryHeap( comparator, events_queue );
    #println( typeof( heap ) );
    push!( heap, ID_START );

    #println( "INIT" );
    fr = FEVERContext( P, Q, eid, vals, events_queue, comparator, heap,
                         prev,
                         cardin( P ), cardin( Q ), false, zero( T )  )
    #println( "INIT" );
    return  fr;
end


function  fever_is_schedule_event( c::FEVERContext{N,T}, id::Int64,
                                   i::Int64, i_is_vert::Bool,
                                   j::Int64, j_is_vert::Bool
                                   )::Bool where {N,T}
    ( ( id > length( c.prev ) ) ||  ( c.prev[ id ] != 0 ) )  &&  return  false

    n_p = c.n_p;
    n_q = c.n_q;

    if  ( ( i  > n_p )  ||  ( j > n_q ) )
        return  false;
    end
    if  ( ( i >= n_p )  &&  ( ! i_is_vert ) )
        return  false
    end
    if  ( j >= n_q )  &&  ( ! j_is_vert )
        return  false
    end
    return  true
end


function  fever_event_value( c::FEVERContext{N,T}, i::Int64,
    i_is_vert::Bool, j::Int64, j_is_vert::Bool,
    id::Int64
) where {N,T}

    ( c.vals[ id ] >= 0 ) &&  return  vals[ id ];
        
    P = c.P
    Q = c.Q
    if  i_is_vert
        if  j_is_vert
            d = Dist( P[ i ], Q[ j ] );
        else
            d = dist_iseg_nn_point( Q[ j ], Q[ j + 1 ], c.P[ i ] );
        end
        c.vals[ id ] = d;
        return  d;
    end

    if  j_is_vert
        d = dist_iseg_nn_point( P[ i ], P[ i + 1 ], Q[ j ] );
        c.vals[ id ] = d;
        return  d;
    end

    println( "Error: This kind of event is not handled yet..." );
    exit( -1 )
    return  0.0;
end


function  fever_schedule_event( c::FEVERContext{N,T},
                                i::Int64, i_is_vert::Bool,
                                j::Int64, j_is_vert::Bool,
                                 prev_id::Int64,
                               ) where  {N,T}
    id = ID_get( c.eid, i, i_is_vert, j, j_is_vert );
    if  ! fever_is_schedule_event( c, id, i, i_is_vert, j, j_is_vert )
        return
    end
    new_val = fever_event_value( c, i, i_is_vert, j, j_is_vert, id );
    if  ( c.f_upper_bound  &&  ( c.upper_bound < new_val ) )
        return;
    end

    ### Also... Blocks the event from being considered again...
    c.prev[ id ] = prev_id;
    push!( c.heap, id );
end


function  fever_is_final_cell( id::Int64, n_p::Int64, n_q::Int64 )
    return   ( ( EID_i( id ) == ( n_p - 1 ) )
               &&  ( EID_j( id ) == ( n_q - 1 ) ) )
end


function  fever_extract_sol_ids( P::Polygon{N,T}, Q::Polygon{N,T},
                                 prevArr::Vector{Int64}
                                 ) where {N,T}
    out_arr = Vector{Int64}();

    curr = ID_END;
    while  ( curr != ID_START )
        push!( out_arr, curr );
        curr = prevArr[ curr ];
    end
    push!( out_arr, curr );

    reverse!( out_arr );

    return  out_arr;
end

function   ID_status( id::Int64 )::Int64
    return   id & 0x1;
end

function   fever_same_status( arr::Vector{Int64}, curr::Int64,
                                  len::Int64 )

    ( curr == len )  &&  return len;
    id = arr[ curr ];
    status = ID_status( id );

    if  ( ID_status( arr[ curr + 1 ] ) != status )
        return  curr;
    end
    t = curr + 1;
    while  ( t <= len )
        if ( ( ID_status( arr[ t ] ) != status )
             ||  ( arr[ t ] == ID_START )
             ||  ( arr[ t ] == ID_END ) )
            break;
        end
        t = t + 1;
    end
    return  t - 1;
end



function   fever_comp_leash( c::FEVERContext{N,T},
                             P::Polygon{N,T},
                             Q::Polygon{N,T},
                             arr::Vector{Int64}
                             ) where {N,T}
    curr::Int64 = 1;
    len::Int64 = length( arr );
    l_max::T = 0.0;
    l_min::T = 0.0;
    eid::EIDCalc = c.eid;

    fc = FeverCoords( 0, false, 0, false );
    fe = FeverCoords( 0, false, 0, false );

    while  curr <= len
        id_curr::Int64 = arr[ curr ];

        ID_get_fields( eid, id_curr, fc );

        if  ( fc.i_is_vert  &&  fc.j_is_vert )
            l_min = max( l_min, Dist( P[ fc.i ], Q[ fc.j ] ) );
            curr = curr + 1;
            continue;
        end

        if  ( fc.i_is_vert  &&  ( ! fc.j_is_vert ) )
            low = curr;
            hi = fever_same_status( arr, curr, len )

            if ( low == hi )
                curr = curr + 1;
                #new_val = dist_iseg_nn_point( q_a, q_b, P[ fc.i ] );
                #@assert( floating_equal( c.vals[ id_curr ], new_val ) );
                l_min = max( l_min, c.vals[ id_curr ] );
                continue;
            end

            ( low > hi )  && @assert( false );

            ID_get_fields( eid, arr[ hi ], fe );

            if  ( fe.j != fc.j )
                println( low, "...", hi );
                for  i in low:hi
                    println( ID_status( arr[ i ] ) );
                end
                @assert( fc.j == fe.j );
            end

            @assert( fc.i <= fe.i );
            #println( " [[[ ", l_min, "...", l_max, " ]]]", fc.i, "...",  fe.i,
            #"       >>>", fc.j );
            q_a = Q[ fc.j ];
            q_b = Q[ fc.j + 1 ];
            l_min, l_max = max_leash( l_min, l_max, q_a, q_b, P, fc.i,
                fe.i );
            #println( "[[[[ ", l_min, "...", l_max, " ]]]" );
            curr = hi + 1;
            continue;
        end

        if  ( ( ! fc.i_is_vert )  &&  (  fc.j_is_vert ) )
            low = curr;
            hi = fever_same_status( arr, curr, len )
            p_a = P[ fc.i ];
            p_b = P[ fc.i + 1 ];
            if ( low == hi )
                curr = curr + 1;
                #new_val = dist_iseg_nn_point( p_a, p_b, Q[fc.j] );
                l_min = max( l_min, c.vals[ id_curr ] );
                #@assert( floating_equal( c.vals[ id_curr ], new_val ) );
                continue;
            end

            ( low > hi ) &&  @assert( false );

            ID_get_fields( eid, arr[ hi ], fe );
            if  ( fe.i != fc.i )
                println( low, "...", hi );
                for  i in low:hi
                    println( ID_status( arr[ i ] ) );
                end
                @assert( fc.i == fe.i );
            end


            @assert( fc.j <= fe.j );
            l_min, l_max = max_leash( l_min, l_max, p_a, p_b, Q, fc.j, fe.j );
            curr = hi + 1;
            continue;
        end

        @assert( false );
    end

    l_max = max( l_min, l_max );
    return  l_min, l_max;
end


function   FEVER_compute_range( P::Polygon{N,T},
                                Q::Polygon{N,T},
                                upper_bound::T
                               ) where {N,T}
    #println( "\n\n\n\n\n" );
    f_debug::Bool = false;
    c::FEVERContext{N,T} = FEVER_Context( P, Q )
    fc = FeverCoords( 0, false, 0, false );

    ## DEBUG code
    #ID_tester( cardin( P ), cardin( Q ) );

    c.f_upper_bound = true;
    c.upper_bound = upper_bound;

    n_pm::Int64 = c.n_p - 1;
    n_qm::Int64 = c.n_q - 1;
    heap = c.heap;
    vals = c.vals;
    eid = c.eid;
    iters = 0;
    id::Int64 = 0;
    while  ! isempty( heap )
        iters = iters + 1;
        #println( "pop..." );
        id = pop!( heap );
        val = vals[ id ];

        if  id == ID_START
            fever_schedule_event( c, 1, false, 1, true, id );
            fever_schedule_event( c, 1, true , 1, false, id );
            continue;
        end

        ID_get_fields( eid, id, fc );

        #println( "K", Int64( fc.i_is_vert ),
        #    Int64(fc.j_is_vert ),
        #    " (", fc.i, ", ", fc.j, ") ", val );

        if  ( fc.i >= n_pm )  &&  ( fc.j >= n_qm )
            # Is it the *final* event?
            if  ( fc.i == c.n_p )  &&  ( fc.j == c.n_q )
                break;
            end
            if  ( ( fc.i == n_pm )  &&  ( fc.j == n_qm ) )
                c.prev[ ID_END ] = id;
                push!( c.heap, ID_END );
                continue;
            end
        end
        fever_schedule_event( c, fc.i+1, true,   fc.j, false, id );

        fever_schedule_event( c, fc.i  , false, fc.j+1, true, id );
    end
    #println( "ITERS: ", iters );

    #println( "extract..." );
    out_arr = fever_extract_sol_ids( P, Q, c.prev );

    #println( "Compute leash..." );
    l_min, l_max = fever_comp_leash( c, P, Q, out_arr )

    #println( "\n\n\n\n" );
    return  l_min,l_max
end


##########################################################

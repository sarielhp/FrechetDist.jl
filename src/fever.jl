##################################################################3
# fever.jl
#
# An attempt to implement a faster version of the VE Retractable
# Fréchet distance. The idea is to use static arrays to minimize
# dynamic algorithms.
##################################################################3


ID_START = 1
ID_END = 2

struct EIDCalc
    nP::Int64
    nQ::Int64
end

function   ID_init( nP::Int64, nQ::Int64 )
    e = EIDCalc( nP, nQ );

    return  e;
end

function   ID_val_inner( i, j, r, nP )
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

    r = ( i_is_vert ) ? 0 : 1 ;
    # i ∈ 1...nP-1, j ∈ 1.. nQ-1

    return  ID_val_inner( i, j, r, e.nP );
end

function   ID_get_fields( e::EIDCalc, id::Int64 )
    i::Int64;

    if  ( id == ID_START)
        return  1, true, 1, true;
    elseif  ( id == ID_END)
        return e.nP, true, e.nQ, true;
    end
    i_is_vert::Bool = ( ( id & 1 ) != 0 )

    ix = id >> 1;
    j = floor( Int64, (ix - 1)  / nP );
    i = ix - j * nP;

    return  i, i_is_vert, j, ! i_is_vert;
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
            push!( v, ID_get( e, i, false, j, true ) );
            push!( v, ID_get( e, i, true, j, false ) );
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
    println( "nP: ", nP, " nQ: ", nQ, "LB: ", nP * nQ * 2, " vs. ", length( v ), "  mx: ", mx );
end

mutable struct  EventsOrder{T} <: Base.Order.Ordering
    values::Vector{T}
end

import Base.Order.lt
function  lt(o::EventsOrder, a, b)
    return  isless( o.values[ a ], o.values[ b ] );
end

mutable struct FEVERContext{N,T}
    P::Polygon{N,T};
    Q::Polygon{N,T};

    eid::EIDCalc;
    vals::Vector{T}   # An array containing all the event values.
    events_queue::Vector{Int64}  # an array containing the IDs of all
                                # the events in the heap
    comparator::EventsOrder
    heap::BinaryHeap{EventsOrder,Int64};

    prev::Vector{Int64}();

    n_p::Int64 #0
    n_q::Int64 #0
    f_upper_bound::Bool
    upper_bound::T
end

function   FEVER_Context(  P::Polygon{N,T},
    Q::Polygon{N,T} ) where {N,T}

    @assert( ( cardin( P ) > 1 )  &&  ( cardin( Q ) > 1 ) );

    eid = ID_init( cardin( P ), cardin( Q ) );
    mx = ID_get_max( eid );
    vals = zeros( T, mx );
    prev = zeros( Int64, mx );

    vals[ ID_START ] = Dist( first( P ), first( Q ) );
    vals[ ID_END ] = Dist( last( P ), last( Q ) );

    comparator = EventsOrder{T}( vals );

    events_queue = Vector{Int64}( ID_START );
    heap = BinaryHeap{EventsOrder,Int64}( comparator, events_queue );

    return FEVERContext( P, Q, eid, val, events_queue, comparator, heap,
                         prev,
                         cardin( P ), cardin( Q ), false, zero( T )  )
end


function  fever_event_value( c::FEVERContext{N,T}, i::Int64,
    i_is_vert::Int64, j::Int64, j_is_vert::Bool,
    id::Int64
) where {N,T}

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


function  f_r_schedule_event( id::Int64, prev_id::Int64,
                              c::FRContext{N,T} ) where  {N,T}
    if  ! is_schedule_event( c.dict, id, c.n_p, c.n_q )
        return
    end
    new_val = ve_event_value( c, id )
    if  ( c.f_upper_bound  &&  ( c.upper_bound < new_val ) )
        return;
    end
    ### Also... Blocks the event from being considered again...
    c.dict[ id ] = prev_id;
    push!( c.heapAlt, (new_val, id ) );
end




function   fever_compute_range( P::Polygon{N,T},
                                Q::Polygon{N,T},
                                upper_bound::T
                               ) where {N,T}
    f_debug::Bool = false;
    c::FEVERContext{N,T} = FEVER_Context( P, Q )
    c.f_upper_bound = true;
    c.upper_bound = upper_bound;

    n_pm::Int64 = c.n_p - 1;
    n_qm::Int64 = c.n_q - 1;
    heap = c.heap;
    vals = c.vals;
    eid = c.eid;
    while  ! isempty( heap )
        id = pop!( heap );
        val = vals[ id ];

        if  id == START_ID
            fever_schedule_event(  eid, 1, false, 1, true, id, c );
            fever_schedule_event( eid, 1, true, 1, false, id, c );
            continue;
        end

        i, i_is_vert, j, j_is_vert = ID_get_fields( eid, id );


        if  ( i >= n_pm )  &&  ( j >= n_qm )
            # Is it the *final* event?
            if  ( i == c.n_p )  &&  ( j == c.n_q )
                break;
            end
            if  is_final_cell( id, c.n_p, c.n_q )
                c.prev[ ID_END ] = id;
                push!( c.heap, ID_END );
                continue;
            end
        end
        fever_schedule_event( i+1, true, j, false, id, c );
        fever_schedule_event( i, false, j+1, true, id, c );
    end

    out_arr = f_r_extract_solution_ids( P, Q, end_id, start_id, c.dict );

    l_min, l_max = compute_leash_from_arr( P, Q, out_arr )

    return  l_min,l_max
end


















ID_tester( 2, 2 );
ID_tester( 100, 200 );

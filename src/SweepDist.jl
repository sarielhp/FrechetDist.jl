###########################################################################
# SweepDist.jl
#
# Implementation of code to compute the sweeping distance. This is a
# similar entity to the continuous dynamic time warping distance,
# except that we use "L_1" metric instead of "L_2" metric to compute
# the distance (this simplifies the calculations), and of course, we
# use our ve-r framework to compute/approximate it, together with the
# regular refinement for monotonicity, and refinement to get
# convergences.
###########################################################################

DictVerticesType = Dict{Int64, TreeVertex};

@with_kw mutable struct SDContext{N,T}
    P::Polygon{N,T};
    Q::Polygon{N,T};
    vertices::DictVerticesType
    handled::DictHandledType;
    prev_map::DictVERType;
    pq::PriorityQueue{Int64,Float64};
    n_p::Int64 = 0;
    n_q::Int64 = 0;
end

function  SDContext(P::Polygon{N,T}, Q::Polygon{N,T}) where {N,T}
    return SDContext( P, Q,
                        DictVerticesType(),
                        DictHandledType(),
                        DictVERType(),
                        PriorityQueue{Int64,Float64}(),
                        cardin( P ), cardin( Q ) );
end


function  SW_new_event_v( _id::Int64, c::SDContext{N,T},
                            val::Float64 ) where {N,T}
    ev = TreeVertex( _id );
    ev.r = ev.val = 0.0;
    return  ev;
end


function  get_event_points( c::SDContext{N,T}, id::Int64 ) where {N,T}
    P = c.P
    Q = c.Q
    i = EID_i( id );
    j = EID_j( id );
    if  EID_i_is_vert( id )
        if  EID_j_is_vert( id )
            return  P[ i ], Q[ j ];
        end

        q = induced_seg_nn_point( Q[ j ], Q[ j + 1 ], P[ i ] );
        return  P[ i ], q;
    end

    if  EID_j_is_vert( id )
        p = induced_seg_nn_point( P[ i ], P[ i + 1 ], Q[ j ] );
        return  p, Q[ j ];
    end

    println( "Error: This kind of event is not handled yet..." );
    return  P[1], Q[1];
end



function  SW_get_event_value( _id::Int64, _prev_id::Int64, prev_val::Float64,
                                c::SDContext{N,T} ) where {N,T}
    ev = TreeVertex( _id );

    o_p, o_q = get_event_points( c, _prev_id );
    n_p, n_q = get_event_points( c, _id );

    return   prev_val + segs_match_price( o_p, n_p, o_q, n_q );
end


function  SW_new_event( id::Int64, id_prev::Int64, value::Float64,
                          c::SDContext{N,T} ) where {N,T}
    ev = TreeVertex( id );
    ev.id_prev = id_prev;
    ev.r = ev.val = value;

    c.vertices[ id ] = ev;

    return  ev;
end

function  SW_schedule_event( id::Int64, id_prev::Int64,
                               c::SDContext{N,T} ) where  {N,T}
    if  ! is_schedule_event( c.handled, id, c.n_p, c.n_q )
        return
    end
    ev_prev::TreeVertex = c.vertices[ id_prev ];
    new_val = SW_get_event_value( id, id_prev, ev_prev.val, c );
    local ev;
    if  haskey( c.vertices, id )
        ev = c.vertices[ id ];
        if  ( ev.val < new_val )  # already better event scheduled.
            return;
        end
        ev.val = new_val;
    else
        ev = SW_new_event( id, id_prev, new_val, c );
    end

    ev.id_prev = id_prev;
    if  haskey( c.pq, id )
        delete!( c.pq, id );
    end
    enqueue!( c.pq, id, new_val );

    return  ev
end



##########################################################################
# Compute Sweeping Distance. Essentially Prim/Dijkstra algorithm
#
# Returns a Morphing that encodes the solution
##########################################################################
function   SweepDist_compute( P::Polygon{N,T}, Q::Polygon{N,T} ) where {N,T}
    f_debug::Bool = false;
    c::SDContext{N,T} = SDContext( P, Q )

    start_id = EID( 1, true, 1, true );

    id_end = EID( c.n_p, true, c.n_q, true );
    start_event = SW_new_event( start_id, start_id, 0.0, c );

    enqueue!( c.pq, start_id, start_event.val );

    end_event::TreeVertex = start_event;
    iters = 0;
    pq = c.pq;
    while  ! isempty( pq )
        tp = peek( pq );
        id = tp[1];
        value = tp[ 2 ];
        dequeue!( pq );

        if  ( haskey( c.handled, id ) )
            continue;
        end

        ev::TreeVertex = c.vertices[ id ];
        c.handled[ id ] = true;
        c.prev_map[ id ] = ev.id_prev;
        iters = iters + 1;

        if  f_debug  &&  ( (iters % 1000000) == 0 )
            print( "Iters :" );
            print_int_w_commas( iters );
            print( "  ", length( pq ) );
            print( "  " );
            print_int_w_commas( c.n_p * c.n_q * 2 );
            println( "" );
        end

        i = EID_i( id );
        j = EID_j( id );

        if  id == start_id
            SW_schedule_event( EID( 1, false, 1, true ), id, c );
            SW_schedule_event( EID( 1, true, 1, false ), id, c );
            continue;
        end

        # Is it the *final* event?
        if  ( i == c.n_p )  &&  ( j == c.n_q )
            end_event = ev;
            break;
        end

        # Is it on the boundary of the final cell?
        if  is_final_cell( ev.id, c.n_p, c.n_q )
            SW_schedule_event( id_end, ev.id, c );
            continue;
        end
        SW_schedule_event( EID( i+1, true, j, false ), id, c );
        SW_schedule_event( EID( i, false, j+1, true ), id, c );
    end

    pes, qes = f_r_extract_solution( P, Q, id_end, c.prev_map );

    morph::Morphing{N,T} = Morphing_init( P, Q, pes, qes );
    morph.iters = iters;

    return  morph
end


function   SweepDist_compute_refine_mono( poly_a::Polygon{N,T},
                                     poly_b::Polygon{N,T} ) where {N,T}
    ell = Polygon_length( poly_a ) + Polygon_length( poly_b );

    P = poly_a;
    Q = poly_b;

    println( "*** P :", cardin( P ) );
    println( "*** Q :", cardin( Q ) );
    
    local  m_new;
    while  true
        m = SweepDist_compute( P, Q );
        if  ( Morphing_is_monotone( m ) )
            return  m;
        end

        m_new =  Morphing_monotonize( m );
        rate = m_new.monotone_err / ell;

        #println( "Rate: ", rate );
        if  ( rate < 0.01 )
            break;
        end
        poly_a_2 = extract_refined_polygon( P, m.pes, 1 );
        poly_b_2 = extract_refined_polygon( Q, m.qes, 1 );

        P = poly_a_2;
        Q = poly_b_2;
    end

    println( "Returning from SweepDist_compute_refine_mono" );
    return  m_new;
end

function   SweepDist_compute_split(
    poly_a::Polygon{N,T},
    poly_b::Polygon{N,T},
    out::Vector{Morphing{N,T}},
    iters::Int64,
    max_vertices::Int64 = 1000000 ) where {N,T}

    m = SweepDist_compute_refine_mono( poly_a, poly_b );
    push!( out, m );
    P = m.P;
    Q = m.Q;
    for  i in 1:iters
        P = Polygon_split_edges( P );
        Q = Polygon_split_edges( Q );
        mi = SweepDist_compute( P, Q );
        #println( "i : ", i );
        push!( out, mi );
        if  ( ( cardin( P ) + cardin( Q ) ) > max_vertices )
            break;
        end
    end
end


###################################################################

"""
    SweepDist_lb_compute

Compute a lower bound on the Sweep distance between the two curves P
and Q. Current implementation is silly, and very similar to the
discrete Frechet calculation.

"""
function   SweepDist_lb_compute( P::Polygon{D,T}, Q::Polygon{D,T}
                           ) where {D,T}

    pq = PriorityQueue{Tuple{Int64, Int64},Float64}();

    n_p = cardin( P );
    n_q = cardin( Q );

    dp = zeros( n_p, n_q );
    dp_dec_i = falses( n_p, n_q );
    handled = falses( n_p, n_q );

    """
        edge_lb

    The edge is (i,j) -> (i+1,j). The function returns the minimum
    elvation fucntion, in the two adjacent cells:

    min( dist( PA[i, i+1], PB[ j-1, j ] ),
         dist( PA[i, i+1], PB[ j, j+1 ] ),
    """
    function edge_lb( PA::Polygon{D,T}, i, PB::Polygon{D,T}, j ) where {D,T}
        #println( "i: ", i, "  [", cardin( PA ), "]" );
        #println( "j: ", j, "  [", cardin( PB ), "]" );
        #Dist( PA[ i ], PB[ j ] );
        lb::Float64 = typemax(Float64);
        if  i == cardin( PA )
            return  0.0;
        end
        #println( "PA[i,i+1]:", PA[i], PA[ i + 1 ]);
        if  j > 1
            d = iseg_iseg_dist( PA[i ], PA[i + 1], PB[ j - 1 ], PB[ j ] );
            lb = min( lb, d );
            #println( "PB[j-1,j]:", PB[j-1], PB[ j ] );
            #println( "d: ", d );
        end
        if  j < cardin( PB )
            d = iseg_iseg_dist( PA[ i ], PA[i + 1], PB[ j ], PB[ j + 1 ] )
            lb = min( lb, d );
            #println( "PB[j-1,j]:", PB[j], PB[ j + 1 ] );
            #println( "d: ", d );
        end
        #println( "LB: ", lb );
        #println( "\n\n\n" );
        return  lb;
    end

    function edge_v( PA::Polygon{D,T}, i::Int64, PB::Polygon{D,T},
                     j::Int64 )::Float64 where {D,T}
        if  ( i == cardin( PA ) )
            return 0.0;
        end
        #println( "edge_lb( ?, ", i, ", ?, ", j, ")" );
        lb = edge_lb( PA, i, PB, j );
        len = Dist( PA[ i ], PA[ i + 1 ] );
        #println( "Len: ", len );
        return   len * lb;
    end


    function  push_value( i::Int64, j::Int64,
        ia::Int64, ja, new_val::Float64 )
        #new_val = val + Dist( P[ ia ], Q[ ja ] )

        if  ( ( dp[ ia, ja ] > 0.0 )
              &&  ( dp[ ia, ja ] <= new_val ) )
            return;
        end
        #println( "(", i, ", ", j, ") -> (", ia, ", ", ja, "): ",
        #    new_val );
        dp[ ia, ja ] = new_val;
        pq[ ia, ja ] = new_val;  # Push to queue
        dp_dec_i[ ia, ja ] = ( i < ia );
    end

    dp[ 1, 1 ] = 0.0;
    enqueue!( pq, (1,1), dp[ 1, 1 ] );

    iters::Int64 = 0;
    while  ! isempty( pq )
        iters = iters + 1;

        ele = peek( pq );
        i::Int64, j::Int64 = ele[ 1 ];
        value = ele[ 2 ];

        #println( pq );
        #println( "LOC: (", i, ": ", n_p, ", ",j ,":", n_q, " ): ",
        #    value );
        #println( "  VALUE   : ", value );
        #println( "  dp[", i, ", ", j, "] : ", dp[ i, j ] );
        #println( pq );

        dequeue!( pq );
        if  handled[ i, j ]
            continue;
        end


        handled[ i, j ] = true;
        if  ( ( i == n_p )  &&  ( j == n_q ) )
            break;
        end


        # Now we need to schedule the next two adjacent entries..
        if  ( i < n_p )
            push_value( i, j, i + 1, j, value + edge_v( P, i, Q, j ) );
        end
        if  ( j < n_q )
            push_value( i, j, i, j+1, value + edge_v( Q, j, P, i ) );
        end
    end

    m = d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
    m.iters = iters;

    m.sol_value = dp[ n_p, n_q ];

    #println( "sol_value: ", m.sol_value );
    #exit( -1 );
    return  m;
end

# SweepDist.jl end of file
#############################################################################

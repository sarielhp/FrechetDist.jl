# Originally contributed by S. Har-Peled
# under MIT License

#-------------------------------------------------------------
# Frechet
#
# A library to compute various versions of the Frechet distance.
#-------------------------------------------------------------

DictVERType = Dict{Int64, Int64};
#DictHandledType = Dict{Int64, Bool};


"""
    EID

A function that encodes a matching "edge" betwee two polygons.
"""
function  EID(
    i::Int64 = 0,
    i_is_vert::Bool = false,
    j::Int64 = 0,
    j_is_vert::Bool = false
)::Int64
    #---
    val::Int64 = ( ( i << 32 )  |  ( j << 2 )
                   |  ( convert( Int64, i_is_vert ) << 1)
                   |  ( convert( Int64, j_is_vert ) )
    );
    return  val;
end

@inline function  EID_i( id::Int64 )::Int64
    return  ( id >> 32 );
end

@inline function  EID_j( id::Int64 )::Int64
    return   ( id & 0xfffffffc ) >> 2;
end

@inline function   EID_ij_status( id::Int64 )::Int64
    return   id & 0x3;
end

@inline function   EID_i_is_vert( id::Int64 )::Bool
    return   ( (id & 0x2) != 0 );
end

@inline function  EID_j_is_vert( id::Int64 )::Bool
    return  ( (id & 0x1) != 0 );
end

#DictVERType = Dict{Int64, TreeVertex};
#HeapVERType = BinaryMinHeap{TreeVertex};

HeapVerticesType = BinaryMinHeap{Tuple{Float64,Int64},};


#@with_kw
mutable struct FRContext{N,T}
    P::Polygon{N,T};
    Q::Polygon{N,T};
    p_offs::Vector{Float64};
    q_offs::Vector{Float64};
#    handled::DictHandledType;
    dict::DictVERType;
    heapAlt::HeapVerticesType;
    f_offsets::Bool ;
    n_p::Int64 #0
    n_q::Int64 #0
    f_upper_bound::Bool
    upper_bound::T
end

function  FR_Context(P::Polygon{N,T}, Q::Polygon{N,T}) where {N,T}
#    d_h = DictHandledType();
    d_ver = DictVERType();
    h = HeapVerticesType();
    return FRContext( P, Q, Vector{Float64}(),  Vector{Float64}(),
                      d_ver,
                      h,
                      false, Int64(cardin( P) ), Int64(cardin( Q )),
                      false, zero(T)  # Upper bound
                      );
end

function  ve_event_value( c::FRContext{N,T}, id::Int64 ) where {N,T}
    P = c.P
    Q = c.Q
    i = EID_i( id );
    j = EID_j( id );
    if  EID_i_is_vert( id )
        if  EID_j_is_vert( id )
            d = Dist( P[ i ], Q[ j ] );
            # Old code: Bug? Missing if check
            if  c.f_offsets
                return  d -  c.p_offs[ i ] - c.q_offs[ j ];
            else
                return  Dist( P[ i ], Q[ j ] );
            end
        end
        d = dist_iseg_nn_point( Q[ j ], Q[ j + 1 ], c.P[ i ] );
        if  c.f_offsets
#            if  ( c.q_offs[ j ] > 0 )
#                println( "Larger than 0!" );
#            end
            #--------------------------------------------------------------
# Old code: Bug?
#            return  d - c.q_offs[ j ];
            return  d - max(c.q_offs[ j ], c.q_offs[ j + 1 ] ) - c.p_offs[ i ];
        else
            return  d;
        end
    end

    if  EID_j_is_vert( id )
        d = dist_iseg_nn_point( P[ i ], P[ i + 1 ], Q[ j ] );
        if  c.f_offsets
# Old code: Bug?
#            return  d - c.p_offs[ i ];
            return  d - max(c.p_offs[ i ], c.p_offs[ i + 1 ] ) - c.q_offs[ j ];
        else
            return  d;
        end
    end

    println( "Error: This kind of event is not handled yet..." );
    exit( -1 )
    return  0.0;
end

#=
function  f_r_new_event( _id::Int64, c::FRContext{N,T} ) where {N,T}
    ev = TreeVertex( _id );
    ev.val = ve_event_value( c, ev.id );
    return  ev;
end
=#

@inline function   is_start_event( id::Int64 )
    return ( ( EID_i( id ) == 1 )
             && ( EID_j( id ) == 1 )
             && EID_i_is_vert( id )
             &&  EID_j_is_vert( id ) )
end

@inline function  is_final_cell( id::Int64, n_p::Int64, n_q::Int64 )
    return   ( ( EID_i( id ) == ( n_p - 1 ) )
               &&  ( EID_j( id ) == ( n_q - 1 ) ) )
end

function f_r_create_event( R::Polygon{N,T}, i::Int64,
                        is_vert::Bool, qr::Point{N,T} ) where {N,T}
    if  ( is_vert )
        ev = EventPoint( R[ i ], i, PT_VERTEX, 0.0 );
        return  ev;
    end

    p::Point{N,T},t::T = iseg_nn_point_ext( R[ i ], R[ i + 1 ], qr );


    #seg = Segment( R[ i ], R[ i + 1 ] );
    #p = nn_point( seg, qr );
    #=
    if  ( isNaN( p ) )
        println( seg );
        println( qr );
        @assert( false );
    end
    =#

    #t::Float64 = Segment_get_convex_coef( seg, p );

    if  ( ! ( 0.0 <= t  &&  t <= 1.0 ) )
        println( R[ i ] );
        println( R[ i + 1 ] );
        println( "qr: ", qr );
        println( "t: ", t );
        @assert( 0.0 <= t  &&  t <= 1.0 );
    end
    if  ( 0 < t < 0.000001 )
        t = 0;
        p = R[ i ];
    end
    if  ( 0.999999 < t < 1.0 )
        t = 1.0;
        p = R[ i + 1 ];
    end
    #        println( "TTT= ", t );

    @assert( ! isNaN( p ) );
    return  EventPoint( p, i, PT_ON_EDGE, t );
end


@inline function  is_schedule_event( dict::DictVERType, id::Int64, n_p::Int64,
                             n_q::Int64 )::Bool
    if  haskey( dict, id )
        return  false
    end
    if  ( EID_i( id )  < n_p )  &&  ( EID_j( id ) < n_q )
        return  true;
    end
    if  ( EID_i( id )  > n_p )  ||  ( EID_j( id ) > n_q )
        return  false;
    end
    if  ( EID_i(id) >= n_p )  &&  ( ! EID_i_is_vert( id ) )
        return  false
    end
    if  ( EID_j( id ) >= n_q )  &&  ( ! EID_j_is_vert( id ) )
        return  false
    end
    return  true
end

@inline function  f_r_schedule_event( id::Int64, prev_id::Int64,
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


function  f_r_extract_solution( P::Polygon{N,T}, Q,
                                end_event_id::Int64,
                                start_id::Int64,
                                dict::DictVERType
                                ) where {N,T}
    #############################################################3
    # Extracting the solution
    #
    # The final event is always the same - final vertex of P versus
    # final vertex of Q.
    pes = Vector{EventPoint{N,T}}();
    qes = Vector{EventPoint{N,T}}();

    idx = end_event_id;
    pex::EventPoint = f_r_create_event( P, EID_i( idx ),
                                        EID_i_is_vert( idx ),
                                        Q[ EID_j( idx) ] );
    qex::EventPoint = f_r_create_event( Q, EID_j( idx),
                                        EID_j_is_vert( idx ),
                                        P[ EID_i( idx ) ] );

    push!( pes, pex );
    push!( qes, qex );

    curr = end_event_id;
    while  ! is_start_event( curr )
        prev = curr;
        curr = dict[ prev ];

#        id = curr.id;
        id_i = EID_i( curr );
        id_j = EID_j( curr );
        pe::EventPoint = f_r_create_event( P, id_i, EID_i_is_vert( curr ),
                                                 Q[ id_j ] );
        qe::EventPoint = f_r_create_event( Q, id_j, EID_j_is_vert( curr ),
                                              P[ id_i ] );
        push!( pes, pe );
        push!( qes, qe );
    end

    reverse!( pes );
    reverse!( qes );

    return  pes, qes;
end

function   eid_same_status( arr::Vector{Int64}, curr::Int64, len::Int64 )
    eid = arr[ curr ];
    status = EID_ij_status( arr[ curr ] );
    t = curr + 1;
    while  ( t <= len )
        if EID_ij_status( arr[ t ] ) != status
            break;
        end
        t = t + 1;
    end
    return  t - 1;
end


function  f_r_extract_solution_ids( P::Polygon{N,T}, Q::Polygon{N,T},
                                    end_event_id::Int64,
                                    start_event_id::Int64,
                                    dict
                                    ) where {N,T}
    out_arr = Vector{Int64}();

    push!( out_arr, end_event_id );

    curr = end_event_id;
    while  ( curr != start_event_id )
        push!( out_arr, curr );
        prev = curr;
        if  ! haskey( dict, prev )
            println( "start_event_id: ", start_event_id );
            println( "end_event_id: ", end_event_id );
            println( "|P|: ", cardin( P ) );
            println( "|Q|: ", cardin( Q ) );
            println( P );
            println( "----------------------------------" );
            println( Q );
            println( "----------------------------------" );
        end
        curr = dict[ prev ];
    end
    push!( out_arr, curr );

    reverse!( out_arr );

    return  out_arr;
end

"""
    max_leash

    Gets a segment on one polygon s_a s_b, and a chain P[low:hi]. Updates the
    min/max estimates from the leash length. Lower bound is simply the distance
    to nearest point on th esegment. The max is the result of brute force
    monotonization.
"""
@inline function    max_leash( l_min::T, l_max::T,
                               s_a::Point{N,T}, s_b::Point{N,T},
                               P::Polygon{N,T}, low::Int64,
                               hi::Int64 ) where {N,T}
    #seg = Segment( s_a, s_b );
    len_seg_sq::Float64 = DistSq( s_a, s_b );
    if  ( len_seg_sq == 0.0 )
        for  j in low:hi
            l_min = max( l_min, Dist( s_a, P[ j ] ) );
        end
        return  l_min, max( l_min, l_max );
    end

    max_t = 0.0;
    max_s = s_a;

    seg_len_sq = DistSq( s_a, s_b )  # Squared distance between s_a and s_b
    s_vec =  sub(s_b, s_a)

    for  j in low:hi
        p = P[ j ];

        sq, new_t = iseg_nn_point_ext_ext( s_a, s_b, p, seg_len_sq, s_vec );

        dst = Dist( sq, p );
        l_min = max( l_min, dst );
        if  ( new_t >= max_t )
            max_t = new_t;
            max_s = sq; #convex_comb( p_a, p_b, max_t );
        end

        l_max = max( l_max, Dist( max_s, p ) );
    end

    return  l_min, l_max;
end

function   compute_leash_from_arr( P::Polygon{N,T},
                                   Q::Polygon{N,T},
                                   arr::Vector{Int64}
                                   ) where {N,T}
    curr::Int64 = 1;
    len::Int64 = length( arr );
    l_max::T = 0.0;
    l_min::T = 0.0;

    while  curr <= len
        eid = arr[ curr ];

        i = EID_i( eid );
        j = EID_j( eid );

        if  EID_i_is_vert( eid )  &&  EID_j_is_vert( eid )
            l_min = max( l_min, Dist( P[ i ], Q[ j ] ) );
            curr = curr + 1;
            continue;
        end

        if  ( EID_i_is_vert( eid )  &&  ( ! EID_j_is_vert( eid ) ) )
            low = curr;
            hi = eid_same_status( arr, curr, len )

            if  low > hi
                curr = curr + 1;
                continue;
            end
            eid_start = arr[ low ];
            eid_end = arr[ hi ];
            q_a = Q[ EID_j( eid_start ) ];
            q_b = Q[ EID_j( eid_start ) + 1 ];

            if ( low == hi )
                curr = curr + 1;
                l_min = max( l_min, dist_iseg_nn_point( q_a, q_b,
                                                P[ EID_i( eid_start ) ] )
                             );
                continue;
            end

            #println( " ((( ", l_min, "...", l_max, " )))", EID_i( eid_start ),
            #    "...",  EID_i( eid_end ),  "    >>>", EID_j( eid_start ) );
            l_min, l_max = max_leash( l_min, l_max, q_a, q_b, P,
                                      EID_i( eid_start ),
                                      EID_i( eid_end ) );
            #println( "(((( ", l_min, "...", l_max, " )))");
            curr = hi + 1;
            continue;
        end

        if  ( ( ! EID_i_is_vert( eid ) )  &&  (  EID_j_is_vert( eid ) ) )
            low = curr;
            hi = eid_same_status( arr, curr, len )

            if  low > hi
                curr = curr + 1;
                continue;
            end
            eid_start = arr[ low ];
            eid_end = arr[ hi ];
            p_a = P[ EID_i( eid_start ) ];
            p_b = P[ EID_i( eid_start ) + 1 ];

            if ( low == hi )
                curr = curr + 1;
                l_min = max( l_min,
                             dist_iseg_nn_point( p_a, p_b,
                                                Q[ EID_j( eid_start ) ] ) );
                continue;
            end

            l_min, l_max = max_leash( l_min, l_max, p_a, p_b, Q, EID_j( eid_start ),
                               EID_j( eid_end ) );
            curr = hi + 1;
            continue;
        end

        @assert( false );
    end

    l_max = max( l_min, l_max );
    return  l_min, l_max;
end


"""
    frechet_ve_r_compute_range

    Return a range of distances. The lower end is the VE-Frechet
    distance, the upper end is a constructive upper bound on the
    Frechet distnace.

"""
function   frechet_ve_r_compute_range( P::Polygon{N,T},
                                       Q::Polygon{N,T},
                                       upper_bound::T
                                       ) where {N,T}
    f_debug::Bool = false;
    c::FRContext{N,T} = FR_Context( P, Q )
    c.f_upper_bound = true;
    c.upper_bound = upper_bound;

    n_pm::Int64 = c.n_p - 1;
    n_qm::Int64 = c.n_q - 1;

    start_id = EID( 1, true, 1, true );

    end_id = EID( c.n_p, true, c.n_q, true );
    f_r_schedule_event( start_id, start_id, c );

    heapAlt = c.heapAlt;
    while  ! isempty( heapAlt )
        tp = pop!( heapAlt );
        id = last( tp );
        value = first( tp );

        i = EID_i( id );
        j = EID_j( id );

        #println( "J", Int64(EID_i_is_vert( id )),
        #     Int64(EID_j_is_vert( id )),
        #    " (", i, ", ", j, ") ", value );

        if  id == start_id
            f_r_schedule_event( EID( 1, false, 1, true ), id, c );
            f_r_schedule_event( EID( 1, true, 1, false ), id, c );
            continue;
        end

        if  ( i >= n_pm )  &&  ( j >= n_qm )
            # Is it the *final* event?
            if  ( i == c.n_p )  &&  ( j == c.n_q )
                break;
            end
            if  is_final_cell( id, c.n_p, c.n_q )
                f_r_schedule_event( end_id, id, c );
                continue;
            end
        end
        f_r_schedule_event( EID( i+1, true, j, false ), id, c );
        f_r_schedule_event( EID( i, false, j+1, true ), id, c );
    end

    out_arr = f_r_extract_solution_ids( P, Q, end_id, start_id, c.dict );

    l_min, l_max = compute_leash_from_arr( P, Q, out_arr )

    return  l_min,l_max
end



##########################################################################
##########################################################################
# Compute retractable Frechet. Essentially Prim/Dijkstra algorithm
#
# Returns a Morphing that encodes the solution
##########################################################################
function   frechet_ve_r_compute_ext( P::Polygon{N,T},
                                     Q::Polygon{N,T},
                                     p_offs::Vector{Float64},
                                     q_offs::Vector{Float64},
                                     f_use_offsets::Bool = false,
                                     ) where {N,T}
    f_debug::Bool = false;
    c = FR_Context( P, Q )

    if  f_use_offsets
        c.p_offs = p_offs;
        c.q_offs = q_offs;
        c.f_offsets = true;
    end

    start_id = EID( 1, true, 1, true );

    end_id = EID( c.n_p, true, c.n_q, true );
    f_r_schedule_event( start_id, start_id, c );

    #end_event_i::TreeVertex = start_event;
    iters = 0;
    heapAlt = c.heapAlt;
    while  ! isempty( heapAlt )
        tp = pop!( heapAlt );
        id = last( tp );
        value = first( tp );
        #if  ( haskey( c.handled, id ) )
        #  continue;
        #end
        #c.handled[ id ] = true;
        iters = iters + 1;

        if  f_debug  &&  ( (iters % 10000) == 0 )
            print( "Iters :" );
            print_int_w_commas( iters );
            print( "  ", length( heap ) );
            print( "  " );
            print_int_w_commas( c.n_p * c.n_q * 2 );
            print( "   dict size: " );
            print_int_w_commas( length( c.dict ) );
            println( "  " );
        end

        i = EID_i( id );
        j = EID_j( id );

        if  id == start_id #is_start_event( id )
            f_r_schedule_event( EID( 1, false, 1, true ), id, c );
            f_r_schedule_event( EID( 1, true, 1, false ), id, c );
            continue;
        end

        # Is it the *final* event?
        if  ( i == c.n_p )  &&  ( j == c.n_q )
            break;
        end

        # Is it on the boundary of the final cell?
        if  is_final_cell( id, c.n_p, c.n_q )
            f_r_schedule_event( end_id, id, c );
            continue;
        end
        f_r_schedule_event( EID( i+1, true, j, false ), id, c );
        f_r_schedule_event( EID( i, false, j+1, true ), id, c );
    end

    #@time
    pes, qes = f_r_extract_solution( P, Q, end_id, start_id, c.dict );

    morph::Morphing{N,T} = Morphing_init( P, Q, pes, qes );
    morph.iters = iters;

    return  morph
end


function   frechet_ve_r_compute( P::Polygon{N,T}, Q::Polygon{N,T} ) where {N,T}
    return frechet_ve_r_compute_ext( P, Q, Vector{Float64}(),
                                     Vector{Float64}(), false );
end




"""
    frechet_width_approx

 2-approximation to the Frechet distance between
   seg = P[first(rng)]-P[last(rng)] and he polygon
   P[rng]
 Here, rng is a range i:j

This function implements a greedy matching alnog the segment - you
move to the cloest point on seg to the current point, making sure
never moving back.

"""
function  frechet_width_approx( P::Polygon{N,T},
                                rng::UnitRange{Int64} = 0:0
                                )::T where {N,T}
    card = cardin( P );
    if  ( card <= 2 )
        return  0;
    end

    if  ( rng == 0:0 )
        rng = 1:card
    end

    if  ( length( rng ) <= 2 )
        return  0;
    end

    s_p = P[ first( rng ) ];
    s_q = P[ last( rng ) ];

    seg_len_sq = DistSq( s_p, s_q )
    s_vec =  sub(s_q, s_p)

    t::Float64 = 0;
    curr::Point{N,T} = s_p;
    leash::Float64 = 0;
    for  i  in  first(rng)+1:last(rng)-1
        q,new_t = iseg_nn_point_ext_ext( s_p, s_q, P[ i ], seg_len_sq, s_vec );
        #q = nn_point( seg, P[ i ] );
        #new_t = Segment_get_convex_coef( seg, q );
        if  ( new_t > t )
            t = new_t;
            curr = q;
        end
        leash = max( leash, DistSq( curr, P[ i ] ) );
    end
    return  sqrt( leash );
end

eachindexButLast(x) =  firstindex(x):lastindex(x)-1

#####################################################################
"""
    frechet_offests

Given P, and subset of indices of P defining a subpolygon, computes a
-approximation to the Frechet distance between every edge of the
subpolygon they define, and corresponding portion of P.  A fast, and
probably decent approximation to the optimal morhping. Finally, set
for every vertex of the subcurve, the frechet distance of this
subcurve. Thus, providing a fast estimates of the offsets.

"""
function    frechet_offsets( P::Polygon{N,T}, p_indices::Vector{Int64}
                             ) where {N,T}
    offs = Vector{Float64}();

    for  i in eachindexButLast( p_indices )
        d = frechet_width_approx( P, p_indices[ i ]:p_indices[ i+1 ] );
        push!( offs, d );;
    end
    return  offs;
end


#####################################################################
"""
    discrete_frechet_sample

Computes the discrete Frechet distance between the two curves.by
sampling them

It first sapmles the two curves rougly uniformly. n is supposed to
be the nubmer of vertices, by the optput might be a sample that is
slightly bigger, as the code tries to ensure the vertices are being
picked.

# Arguments

- `n`: number of vertices to add when "refining" the two curves. The
  number of vertices computed might be larger (but hopefully not much
  larger).

- f_lopt = true: Use standard discrete Frechet or the retractable
  version.

"""
#####################################################################
function   discrete_frechet_sample( polya::Polygon{N,T},
                                    polyb::Polygon{N,T},
                                    n::Int64,
                                    f_lopt::Bool = true
                                    )    where {N,T}
#    polya,polyb = example_4()
    lena = polygon.total_length( polya );
    lenb = polygon.total_length( polyb );

#    n::Int64 = 100;
    delta = (lena+ lenb)/ n;

    na::Int64 = ceil( lena / delta ) +2;
    nb::Int64 = ceil( lenb / delta ) +2;

    qa = Polygon_sample_uniformly( polya, na );
    qb = Polygon_sample_uniformly( polyb, nb );

    if  ( f_lopt )
        morph = frechet_d_r_compute( qa, qb );
    else
        morph = frechet_d_compute( qa, qb );
    end
#    output_morphing( morph, filename, polyb );

    return  morph
end


"""
    frechet_c_mono_approx_subcurve

Approximates the Frechet distance between a curve (P) and subcurbe (Q). Here,
Q vertices are the vertices of P specifieid by p_indices. That is
Q[ i ] =P[ p_indices[ i ]].

"""
function    frechet_c_mono_approx_subcurve(
    P::Polygon{N,T},
    Q::Polygon{N,T},
    p_indices::Vector{Int64} ) where {N,T}

    offsets = zeros( T, length( p_indices ) );
    #println( "#P: ", cardin( P ), "#Q: ", cardin( Q ), "  #qind: ",
    #         length( p_indices ) );

    @assert( cardin(Q) == length( p_indices ) );

    pes = Vector{EventPoint{N,T}}();
    qes = Vector{EventPoint{N,T}}();

    for  i  in  eachindexButLast( p_indices )
        pind = p_indices[ i ];
        pind_next = p_indices[ i + 1 ];
        push!( pes, EventPoint( P[ pind ], pind, PT_VERTEX, 0.0 ) );
        push!( qes, EventPoint( Q[ i ], i, PT_VERTEX, 0.0 ) );

        max_d::Float64 = 0;
        max_t::Float64 = 0;
        seg = Segment( Q[ i ], Q[ i + 1 ] );
        for  j in (pind + 1):(pind_next - 1)
            p = P[ j ];
            q = nn_point( seg, p );
            new_t = Segment_get_convex_coef( seg, q );
            if  ( new_t >= max_t )
                max_t = new_t;
            else # new_t < max_t
                q = segment.at( seg, max_t )
            end
            max_d = max( max_d, Dist( P[ j ], q ) );
            push!( pes, EventPoint( P[ j ], j, PT_VERTEX, 0.0 ) );
            push!( qes, EventPoint( q, i, PT_ON_EDGE, max_t ) );
        end
        offsets[ i ] = max( offsets[ i ], max_d );
        offsets[ i + 1 ] = max( offsets[ i + 1 ], max_d );
    end

    push!( pes, EventPoint( last( P ), cardin( P ), PT_VERTEX, 0.0 ) );
    push!( qes, EventPoint( last( Q ), cardin(Q), PT_VERTEX, 0.0 ) );

    return  Morphing_init( P, Q, pes, qes ), offsets;
end


"""
    add_points_along_seg

# Arguments

seg: t[1], t[2], .... , t[end]

     t[...]: Time stamps for points on the segment (all between 0 and 1)

# Details

Add points along seg at time stamps specified by t, and add points in
the middle between them (i.e., refining segment). These new points are
added the output polygon pout.

Note, that one can set n_mid_points to be larger than 1, but it does
not seem to be necessary/benefitial from the experimentation we did.

Output:
     Sorted( t[1]/2, t[1], (t[1]+t[2])/2, t[2],...)
"""
function   add_points_along_seg( pout::Polygon{N,T},
                                 seg::Segment{N,T},
                                 times::Vector{Float64},
                                 n_mid_points::Int64 = 1 ) where {N,T}
    out_times = Vector{Float64}();

    function  push_time( tm::Float64 )
        if  ( tm <= 0.0 )  ||  ( tm >= 1.0 )
            return
        end

        push!( out_times, tm );
    end

    sort!( times )
    unique!( times )
    times_length = length(times);

    if  isempty( times )
        return
    end

    if  ( n_mid_points > 0 )
        push_time(  times[ 1 ] / 2.0 );
    end

    for  i in 1:length(times)
        t = times[ i ];
        push_time( t );
        if   ( i < times_length )  &&  ( n_mid_points > 0 )
            K = n_mid_points + 1
            #println( "K : ", K );
            #if  ( times[ i + 1 ] > times[ i ] )  continue  end;
            for  pos in 1:(K-1)
                coef::Float64 = pos / K;
                tm = times[ i ] * (1.0 - coef) + times[ i + 1 ] * coef
                push_time( tm );
            end
        end
    end
    if  ( n_mid_points > 0 )
        push_time( ( 1.0 + last(times) ) / 2.0 );
    end

    if  isempty( out_times )
        return;
    end

    sort!( out_times )
    unique!( out_times )

    prev::Float64 = -1.0;
    dlt = 0.09;
    for  tm in out_times
        if  ( ( tm <= 0.0 )  ||  ( tm >= 1.0 ) )
            continue;
        end
        if  ( tm < ( prev + dlt ) )
            continue;
        end
        prev = tm;
        push_smart!( pout, segment.at( seg, tm ) )
    end
end


function  is_monotone_inc( nums::Vector{Float64} )
    len = length( nums );
    if  ( len <= 1 )
        return  true
    end
    for  i in 1:(len-1)
        if  nums[ i ] > nums[ i + 1 ]
            return  false
        end
    end

    return  true
end


################################################################
"""
    extract_refined_polygon( P, pes, num_points_to_add, delta )

    Refine edges along which pes (the sequence of matching points
    along P is not monotone, so that the returned polygon has
    vertices, which hopefully force the Frechet morphing computed in
    the next iteration to be monotone.

# Arguments

  `P`  : the polygo
  `pes`: Sequence of events long the polygon P
  `num_points_to_add`: how many points to add between points of
                       non-monotonicity. The default is 1. Having a
                       larger value speeds up the process of getting
                       rid of non-monotonicity, but the polygons grow
                       faster.
"""
function  extract_refined_polygon( poly::Polygon{N,T},
                                   s::Vector{EventPoint{N,T}},
                                   points_to_add::Int64 = 1,
                                   delta::Float64 = 0.0 ) where {N,T}
    pout = Polygon{N,T}()
    opnts = Vector{EventPoint{N,T}}();

    len = length( s );
    push!( pout, s[1].p );

    i = 2;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
            push_smart!( pout, s[i].p );
            i = i + 1;
            continue;
        end

        j = i
        loc = s[ i ].i;
        while  ( ( j < len )
            &&  ( s[ j + 1 ].type == PT_ON_EDGE )
            &&  ( s[ j + 1 ].i == loc ) )
            j = j + 1
        end

        # s[i:j]: The sequence of points on the same edge
        times = Vector{Float64}();
        push!( times, s[ i ].t );
        for  k in (i+1):j
            if  s[ k ].t != last( times )
                push!( times, s[ k ].t );
            end
        end

        lnx = Dist( poly[ loc ], poly[ loc + 1 ] );

        if  (  is_monotone_inc( times )
               ||  ( lnx < (delta / 2.0 ) ) )
            for  k in i:j
                push_smart!( pout, s[k].p );
            end
            i = j + 1
            continue
        end

        #=
        diff = 0.0;
        min_t = max_t = times[ 1 ];

        for  i in 2:length(times)
            max_t = max( times[ i ], max_t );
            if  ( times[ i ] < max_t )  &&  ( ( max_t - times[ i ] ) > diff )
                min_t = times[ i ];
                diff = max_t - times[ i ];
            end
        end
        =#
        #println( "diff: ", diff );
        #println( "min_t: ", min_t );

        #println( times );
        #exit( -1 );
        n_times::Int64 = points_to_add;
        if  ( delta > 0.0 )
            n_times = round( Int64, lnx / delta );
            n_times = max( min( n_times, 3 ), 1 );
        end
        seg = Segment( poly[ loc ], poly[ loc + 1 ] );
        add_points_along_seg( pout, seg, times, n_times );

        i = j + 1
    end

    return  pout
end


function  frechet_refinement( P::Polygon{N,T}, Q::Polygon{N,T},
                              m::Morphing{N,T} ) where {N,T}
    poly_a_2 = extract_refined_polygon( P, m.pes, 1 );
    poly_b_2 = extract_refined_polygon( Q, m.qes, 1 );

    return  poly_a_2, poly_b_2;
end

function  frechet_mono_via_refinement_ext( Pa::Polygon{N,T}, Qa::Polygon{N,T},
                                           out::Vector{Morphing{N,T}},
                                           f_snapshots::Bool,
                                           approx::Float64 = 1.00001
                                     )  where {N,T}
    @assert( approx > 1.0 );
    f_debug::Bool = false;

    mm::Morphing{N,T} = Morphing_empty( Pa, Qa );
    m::Morphing{N,T} = Morphing_empty( Pa, Qa );
    f_exact::Bool = false;
    f_refinement::Bool = false;
    P = Pa;
    Q = Qa;
    f_debug && println( "fmvr before while..." );
    f_debug && println( "APPROX : ", approx );
    while  true
        f_debug && println( "\n<<<<<<< " );
        f_debug && println( "Before ve_r_compute( ", cardin(P), ", ",
                            cardin( Q ) );
        m = frechet_ve_r_compute( P, Q )
#        f_debug  && println( "After... m.leash : ", m.leash );
        mm = Morphing_monotonize( m );

        if  ( f_snapshots )
            push!( out, m );
        end

        fr_retract = m.leash;
        fr_r_mono = mm.leash;

        if ( f_debug )
            println( "ve_r            : ", m.leash );
            println( "ve_r.mono       : ", fr_r_mono );
            println( "ve_r*approx     : ", approx * fr_retract );
        end

        if  ( fr_r_mono == fr_retract )
            f_exact = true;
            break;
        end
        if  ( fr_r_mono <= approx * fr_retract )
            f_exact = false;
            break;
        end
        f_refinement = true;

        poly_a_2 = extract_refined_polygon( P, m.pes, 1 );
        poly_b_2 = extract_refined_polygon( Q, m.qes, 1 );

        P = poly_a_2;
        Q = poly_b_2;
    end

    # We potentially need to flatten the morphing back to the original
    # curves...
    if   ( ( cardin( mm.P ) != cardin( Pa ) )
           ||  ( cardin( mm.Q ) != cardin( Qa ) ) )
        prm = Morphing_extract_prm( mm );
        mm = Morphing_from_prm( prm, Pa, Qa );
    end

    mm.lower_bound = m.leash;

    #println( "RETURNING!" );
    return  mm, f_exact, P, Q
end


function  frechet_mono_via_refinement_delta( Pa::Polygon{N,T},
    Qa::Polygon{N,T},
    delta::Float64,
    f_fix_prm::Bool = false
)  where {N,T}
    f_debug::Bool = false;

    P = Pa;
    Q = Qa;
    #println( "ROUNDS\n" );
    local  m, mm;
    count::Int64 = 0;
    for  count  in 0:64
        #println( count, "  ", length( P ), "  ", length( Q ) );
        #( count > 0 )  &&  println( "FVER round: ", count );

        #@time m_old = frechet_ve_r_compute( P, Q )
        m = FEVER_compute( P, Q );
        #=
        println( "---" );
        
        if  ( ! fp_equal( m_old.leash, m.leash ) )
            @assert( m_old.leash == m.leash );
        end
        =#
        
        mm = Morphing_monotonize( m );

        res = 3;
        if  ( delta > 0.0 )
            gap = ( mm.leash - m.leash ) / delta;
            if  ( gap > 0.0 )
                res = 1 + (max( 1, round(Int64, (log(gap)/log(2) ) ) ))^2
            end
        else
            gap = delta;
        end

        ( ( mm.leash - m.leash ) <= delta )  &&  break;
        #println( count, ": X GAP: ", gap, "   res: ", res );
        #println( count, ":   res: ", res );

        #times::Int64 = 4;# min( round(Int64,
        #   (mm.leash - m.leash) / delta ), 7 );;
        # res was 3
        poly_a_2 = extract_refined_polygon( P, m.pes, res, delta/11 );
        poly_b_2 = extract_refined_polygon( Q, m.qes, res, delta/11 );

        #        println( "P len: ", length( P ) );
        #println( "P_a_2 len: ", length( poly_a_2 ) );

        P = poly_a_2;
        Q = poly_b_2;
    end

    # We potentially need to flatten the morphing back to the original
    # curves...
    if   f_fix_prm  &&  ( ( cardin( mm.P ) != cardin( Pa ) )
                          ||  ( cardin( mm.Q ) != cardin( Qa ) ) )
        prm = Morphing_extract_prm( mm );
        mm = Morphing_from_prm( prm, Pa, Qa );
    end

    mm.lower_bound = m.leash;

    return  mm, P, Q
end



"""
    frechet_mono_via_refinement( P, Q, approx )

    Computes the "true" monotone Frechet distance between P and Q,
    using the ve_r algorithm. It does refinement, to add vertices if
    needed to get monotonicity. Note, that there is an eps parameter -
    for real inputs you can set it quite small. In any case, in the
    worst case it only approximates the Frechet (monotone)
    distance. It returns a morphing, and a boolean that is true if the
    result is the true Frechet distance.

    Observe, that this function might be too slow if the two curves
    are huge, since it does not simplify them before computing the
    distance.

# Returns

Returns a 4-tuple  ( m, f, P, Q ):
`m`: Morphing realizing result.
`f`: Is distance returned is the exact continuous Frechet distance.
`P`: Refined first input polygon
`Q`: Refined second input polygon

'out' if specified returns the middle morphings used in computing the
      solution.
"""
function  frechet_mono_via_refinement( Pa::Polygon{N,T}, Qa::Polygon{N,T},
                                       approx::Float64 = 1.00001
                                     )  where {N,T}
    ot::Vector{Morphing{N,T}} = Vector{Morphing{N,T}}();
    return  frechet_mono_via_refinement_ext( Pa, Qa, ot, false, approx );
end

###########################################################################
# Compute the morphing src(v) -> trg(u): That is u(v(t))
###########################################################################

"""
    frechet_dist_upper_bound

Returns a rough upper bound on the Frechet distance between the two
curves. This upper bound is on the continuous distance. No guarenteee
on how bad the approximation is. This is used as a starting point for
real approximation of the Frechet distance, and should not be used
otherwise.

"""
function  frechet_dist_upper_bound( poly_a::Polygon{N,T},
                                    poly_b::Polygon{N,T}
                                   ) where {N,T}
    w_a = frechet_width_approx( poly_a );
    w_b = frechet_width_approx( poly_b );
    if  ( ( cardin( poly_a ) < 2 )
          ||  ( cardin( poly_a ) < 2 ) )
        return  w_a + w_b;
    end
    w = max( Dist( poly_a[ 1 ], poly_b[ 1 ] ),
             Dist( last( poly_a ), last( poly_b ) ) )
    return  w_a + w_b + w;
end


"""
    frechet_ve_r_mono_compute

Computing the ver frechet distance, and then monotonize it. No
guarentee as far as the quality of the distance this realizes...

"""
function  frechet_ve_r_mono_compute( poly_a::Polygon{D,T},
                                     poly_b::Polygon{D,T} ) where {D,T}
    m = frechet_ve_r_compute( poly_a, poly_b );
    mm = Morphing_monotonize( m )

    return   mm;
end


"""
    frechet_c_approx

Approximates the continuous Frechet distance between the two input
curves. Returns a monotone morphing realizing it.

# Arguments

- `approx` : The output morhing has Frechet distance <= approx*optimal.

Importantly, approx can be larger than 2, if you want a rough
approximation.

The quality of approximation is available at ret.ratio. Thus,
ret.leash/ret.ratio is a lower bound on the Frechet distance, while
ret.leash is an upper bound.

"""
function  frechet_c_approx( poly_a::Polygon{N,T},
    poly_b::Polygon{N,T}, approx::Float64 )::Morphing{N,T} where {N,T}

    f_debug::Bool = false;

    t_approx::Float64 = approx;
    @assert( ( cardin( poly_a ) > 1 )  &&  ( cardin( poly_b ) > 1 ) )
    @assert( approx > 1.0 )

    f_do_one_round::Bool = true;
    d = frechet_dist_upper_bound( poly_a, poly_b )

    f_debug  &&  println( "DIST :", d );

    # r: Radius of simplification allowed
    r::Float64 = d/( approx + 4.0 ); #max( d/4, d * min( approx, 1.0 ) );

    P::Polygon{N,T} = poly_a;
    Q::Polygon{N,T} = poly_b;

    p_indices = Vector{Int64}();
    q_indices = Vector{Int64}();

    m = Morphing_empty( poly_a, poly_b );
    while  ( true )
        while  ( (r >= ( d / (approx + 4.0) ) )
                 ||  f_do_one_round )
            r = r / 2.0;
            P, p_indices = Polygon_simplify_ext( poly_a, r )
            Q, q_indices = Polygon_simplify_ext( poly_b, r )

            if  ( (2 * cardin( P ) > cardin( poly_a ) )
                  &&  (2 * cardin( Q ) > cardin( poly_b ) ) )
                #println( "DDDD" );
                f_debug  &&  println( "VVV Calling frechet_mono_via_refinement" );
                m = frechet_mono_via_refinement( poly_a, poly_b )[ 1 ];
                m.ratio = 1.0;
                m.lower_bound = m.leash;
                #                exit( -1 ); # DDDD
                return  m;
            end

            @assert( ( cardin( P ) > 1 )  &&  ( cardin( Q ) > 1 ) )

            f_debug &&  println( "cardins: ", cardin( P ), " , ",
                                 cardin( Q )  );
            t_approx = 3.0/4.0 + t_approx / 4.0
            m = frechet_mono_via_refinement( P, Q, t_approx )[1];
            d = m.leash;  #        frechet_ve_r_compute_dist( P, Q )
            f_debug  &&  println( "d : ", d );
            if  ( fp_equal( d, 0 ) )
                break;
            end;
            f_debug  &&  println( "ddd = ", d );
            f_debug  &&  println( "after" );

            # d - 2r: Lower bound on the real Frechet distance.
            if  f_debug
                println( "d: ", d, "  r: ", r, "    #P: ",
                    cardin( P ),
                    " Q: ", cardin( Q ) )
            end
            f_do_one_round = false;
        end

        #mm = Morphing_monotonize( m );
#        println( "mono a" );
        #m_p = frechet_ve_r_mono_compute( poly_a, P )
        #m_p = frechet_ve_r_mono_compute( poly_a, P )
        f_debug  &&  println( "=== f_c_dist( poly_a, P )" );
        m_p = frechet_c_mono_approx_subcurve( poly_a, P, p_indices )[ 1 ];
#        println( "mono b" );

        #m_q = frechet_ve_r_mono_compute( Q, poly_b )
        f_debug  &&  println( "f_c_dist( poly_b, Q )" );
        m_q = frechet_c_mono_approx_subcurve( poly_b, Q, q_indices )[ 1 ];
        Morphing_swap_sides!( m_q );

        #err = max( m_p.leash, m_q.leash );
        #println( "Error: ", err );

        # m = m_p * mm * m_q
        f_debug && println( "Combine m_p, m" );
        mmu = Morphing_monotonize( Morphing_combine( m_p, m ) );
        f_debug && println( "Combine mmu,m_q" );
        mm_out = Morphing_combine( mmu, m_q );

        f_debug && println( "Verification..." );
        f_debug && Morphing_verify_valid( m_p );
        f_debug && Morphing_verify_valid( m );
        f_debug && Morphing_verify_valid( m_q );
        f_debug && Morphing_verify_valid( mmu );

        total_err = m_p.leash + m_q.leash;
        lbx = (m.leash / t_approx) - total_err;
        if  ( d == 0.0 )
            ratio = 1.0;
        else
            ratio = mm_out.leash / lbx;
        end
        f_debug &&  println( "Ratio :", ratio );
        if  ratio <=  approx
            if f_debug
                println( "-----------------------------------" );
                println( "Leash out: ", mm_out.leash );
                println( "Leash lower bound: ", d - 2r );
                println( "approx ratio          : ",ratio  );
                println( "Required approximation: ", 1.0 + approx );
                println( "||P||: ", polygon.total_length( poly_a ),
                    "  ||Q||: ",
                    polygon.total_length( poly_b ) );
            end
            mm_out.ratio = ratio;
            mm_out.lower_bound = m.leash / ratio;
            return  mm_out;

            #            break;
        end
        f_do_one_round = true;
        if  f_debug
            println( "Another final round? ", ratio );
        end
    end
end


function  count_below_zero( arr )
    count = 0;
    for  x in arr
        if  x <= 0.0
            count = count + 1;
        end
    end
    return  count;
end


function  Polygon_simplify_radii_ext( P::Polygon{N,T}, r::Vector{T}
                                      ) where {N,T}
    f_exact::Bool = false;
    #f_debug  &&  println( "#P :", cardin( P ) );
    #f_debug  &&  println( "#r :", length( r ) );

    PS, p_indices = Polygon_simplify_radii( P, r );
    #println( "@!@! |PS|: ", cardin( PS ) );
    @assert( length( p_indices ) == cardin( PS ) );
    if   ( cardin( PS ) <= ( cardin( P ) / 2 ) )
        @assert( cardin( PS ) == length( p_indices ) );
        return  PS, p_indices, false
    end

    p_indices = Vector{Int64}();
    for  i in 1:cardin( P )
        push!( p_indices, i )
    end
    #println( "______________________________" );
    @assert( cardin( P ) == length( p_indices ) );
    return  P, p_indices, true
end


function  propogate_mins( qz, rounds::Int64 = 1 )
    a = deepcopy( qz );
    for  r  in 1:rounds
        for  i  in  eachindex( a )
            if  ( i == 1 )  ||  ( i == length( qz ) )
                continue;
            end

            a[ i ] = min( qz[ i - 1 ], qz[ i ], qz[ i + 1 ] );
        end

        for  i  in  eachindex( a )
            qz[ i ] = a[ i ];
        end
    end
end

function  eq( a::Float64, b::Float64, tolerance::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (tolerance* (abs( a)  + abs(b) ))
end


function  get_refinement_map( PSR::Polygon{N,T}, PS::Polygon{N,T} ) where {N,T}
    @assert( cardin( PSR ) >= cardin( PS ) );

    len = cardin( PSR );
    lenPS = cardin( PS );
    mp = zeros(Int64, len );

    curr = 1
    i = 1
    while  ( curr <= len )
        if  ( i < lenPS )  &&  ( PSR[ curr ] == PS[ i + 1 ]  )
            i = i + 1;
        end
        mp[ curr ] = i;
        curr = curr + 1 ;
    end

    return  mp;
end


function    offsets_copy_map( map_PS, offsets )
    offs = zeros( Float64, length( map_PS ) );
    for  i in eachindex( map_PS )
        offs[ i ] = offsets[ map_PS[ i ] ];
    end
    return  offs;
end


function  simplify_morphing_sensitive( m::Morphing{D,T},
                                       factor::Float64 ) where {D,T}
    pl, ql = Morphing_extract_vertex_radii( m );

    lower_bound = m.lower_bound;

    println( "LLLlower_bound: ", lower_bound );
    pz = ( ( lower_bound * ones( length( pl ) ) ) - pl ) / factor
    qz = ( ( lower_bound * ones( length( ql ) ) ) - ql ) / factor

    propogate_mins( pz, 1 );
    propogate_mins( qz, 1 );

    #println( pz );

    PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext( m.P, pz );
    QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext( m.Q, qz );

    return PS, QS;
end


###########################################################################
"""
    frechet_c_compute

Compute the exact continuous (monotone) Frechet distance between the
two polygons. It should be reasonably fast.

This function is somewhat slower than the approximate versions. Use it
only if you really want the exact answer. Consider using
frechet_continous_approx instead.

# Details

This works by first computing a very rough approximation, followed by
distance senstiave simplification of the curves. It then compute the
monotone fr_ve_r distance between the simplified curves, and it
combine it to get a distance between the two original cuves. It makre
sure the answers are the same, otherwise, it repeates with a finer
simplification/approximation till they are equal.

Finally, the algorithm uses the fr_ve_r_with_offests distance between
the two simplified curves to comptue a lower bound, and make sure this
is equal to the Frechet distance computed. If they are equal, then the
upper/lower bounds on the Frechet distance of the two curves are the
same, which implies that the computed distance is indeed the desired
Frechet distance.

# More details

To really ensure converges, the monotone distance computed between the
simplification is computed using refinement, so tha the ve_r distance

"""
function  frechet_c_compute( P::Polygon{N,T},
                             Q::Polygon{N,T},
                             f_accept_approx::Bool = true
                             )  where {N,T}
    f_debug::Bool = false;

    # The parameters that can be finetunes
    # 2.0, 8.0, 4.0, 10.0 =>  8.74 seconds
    aprx_refinement::Float64 = 2.0;
    factor::Float64          = 4.0
    factor_scale::Float64    = 1.3;
    approx_scale::Float64    = 10.0;

    ##################################################################
    # Tolerance: If |a-b| <= tolerance*( |a| + |b| ) then we consider
    # numbers to be equal. In a perfect world tolerance would be zero.
    # But floating point issues require us to pick some value.
    ##################################################################
    tolerance::Float64  = 0.00001;
    min_approx::Float64 = 1.00000000001;

    f_debug  &&  println( "\n\n\n\n\n\n" );
    f_debug  &&  println( "P#", cardin( P ) )
    f_debug  &&  println( "Q#", cardin( Q ) )

    mf = frechet_c_approx( P, Q, aprx_refinement );
    #println( "=#= P :", cardin( mf.P ) );
    #println( "=#= Q :", cardin( mf.Q ) );
    ratio_2 = mf.ratio;

    len_a = polygon.total_length( P );
    len_b = polygon.total_length( Q );
    f_debug  &&  println( "#", cardin( P ) )
    if  ( mf.leash == 0 )
        return  mf
    end
    approx = min( 1.0 + (len_a + len_b ) / (100.0 * mf.leash ), 1.1 );
    if   f_debug
        println( "" );
        println( "" );
        println( "Ratio         : ", ratio_2 );
        println( "len_a         : ", len_a );
        println( "len_b         : ", len_b );
        println( "2approx_leash : ", mf.leash );
        println( "" );
        println( "" );
        println( "approx rec :", approx );
        println( "" );
        println( "" );
    end

    f_debug && println( approx, " approximation..." );

    if  ratio_2 <= approx
        m = mf;
        #println( "___ #P: ", cardin( m.P ) );
        ratio = ratio_2;
    else
        f_debug  &&  println( "freceht_c_approx( ", approx, ") " );
        m = frechet_c_approx( P, Q, approx );
        #println( "_=_ #P: ", cardin( m.P ) );
        ratio = m.ratio
        f_debug  &&  println( "freceht_c_approx( ", approx, ")...done" );
    end
        #        frechet_ve_r_mono_approx( P, Q, eps );

    if  ( f_debug )
        println( "Approximation computed..." );
        println( "m.leash : ", m.leash );
        println( "" );
    end

    pl, ql = Morphing_extract_vertex_radii( m );
    f_debug  &&  println( "### m.P:", cardin( m.P ) );
    f_debug  &&  println( "### pl :", length( pl ) );
    f_debug  &&  println( "### P :", cardin( P ) )

    lower_bound = m.leash / ratio;
    if  f_debug
        println( "Before while loop..." );
        println( "LOWER BOUND: ", lower_bound );
        println( "upper BOUND: ", m.leash );
    end
    f_debug  &&  println( "\n\n\n\n\n\n" );

    round::Int64 = 0;
    while  true
        round = round + 1;
        aprx_refinement = max( aprx_refinement, min_approx );

        f_debug  &&  println( "-------------------------------------------" );
        #f_debug  &&  println( "A#", cardin( P ) )
        #f_debug  &&  println( "B#", length( pl ) )
        f_debug  &&  println( "Factor: ", factor,
                              "  Approx : ", aprx_refinement );
        pz = ( ( lower_bound * ones( length( pl ) ) ) - pl ) / factor
        qz = ( ( lower_bound * ones( length( ql ) ) ) - ql ) / factor

        propogate_mins( pz, 8 );
        propogate_mins( qz, 8 );


        #propogate_negs( pz );
        #propogate_negs( qz );
        # f_deprintln( pz );

        if  f_debug
            p_count = count_below_zero( pz );
            q_count = count_below_zero( qz );

            println( "pz count <0: ", p_count, " / ", length( pz ) );
            println( "qz count <0: ", q_count, " / ", length( qz ) );
        end
        f_PS_exact::Bool = false;
        f_QS_exact::Bool = false;

        PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext( P, pz );
        QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext( Q, qz );

        #println( "### #QS : ", cardin( QS ), "   q_ind: ",
        #         length( q_indices ) );
        #println( "### #PS : ", cardin( PS ), "   p_ind: ",
        #          length( p_indices ) );

        if  f_debug
            println( "p_count <0: ", p_count, " / ", length( pz ) );
            println( "q_count <0: ", q_count, " / ", length( qz ) );
            println( "|PS| = ", cardin( PS ), "  |P| = ", cardin( P ) );
            println( "|QS| = ", cardin( QS ), "  |Q| = ", cardin( Q ) );
            println( "Computing radii simplified Frechet distance..." );
        end

        # If the PS and QS are really large, this is a waste of
        # time. Just compute the Frechet distance using
        # monotonization.
        if  ( ( 2 * cardin( PS ) > cardin( P ) )
            &&  ( 2 * cardin( QS ) > cardin( Q ) ) )
            f_debug  &&  println( "Simpl. not profitable, doing mono" );
            return   frechet_mono_via_refinement( P, Q )[ 1 ];
        end

        #    m_mid = frechet_ve_r_mono_compute( PS, QS  );
        f_debug  &&  println(  "frechet_mono_via_refinement( PS, QS )" );
        m_mid, _f_exact, PSR, QSR = frechet_mono_via_refinement( PS, QS,
            max( tolerance + 1.0, aprx_refinement ) );

        f_debug  &&  println( "frechet mono via refinement computed" );

        # BUG FIX: Used the mono computation instead of the subcurve code,
        # which is much faster...
        #println( "_#PS: ", cardin( PS ) );
        m_a, offsets_a = frechet_c_mono_approx_subcurve( P, PS, p_indices );

        # OLD CODE
        #   m_a = frechet_ve_r_mono_compute( P, PS );
        mmu = Morphing_combine( m_a, m_mid );
        #f_debug  &&  println( "ve_r_mono( QS -> Q )" );

        # BUG FIX: Used the mono computation instead of the subcurve code,
        # which is much faster...
        #println( "_#QS: ", cardin( QS ) );
        m_b, offsets_b = frechet_c_mono_approx_subcurve( Q, QS, q_indices );
        Morphing_swap_sides!( m_b );

        # OLD CODE
        #  m_b = frechet_ve_r_mono_compute( QS, Q );
        mw = Morphing_combine( mmu, m_b );

        #f_debug  &&  println( "CARDIN(P): ", cardin( mw.P ) );
        #f_debug  &&  println( "CARDIN(Q): ", cardin( mw.Q ) );

        f_debug  &&  println( "eq?  ", m_mid.leash, " = ", mw.leash );

        # is there is NO hope we got the optimal solution?
        if  ( ! eq( m_mid.leash, mw.leash, tolerance ) )
            lb = mw.leash - 2.0*m_a.leash - 2.0*m_b.leash;
            lower_bound = max( lower_bound, lb );
            if  f_debug
                println( "!!! Ratio        : ",
                    fp_ratio( m_mid.leash, mw.leash ) );
                println( "!!! m_mid.leash  : ", m_mid.leash );
                println( "!!! mw.leash     : ", mw.leash );
                println( "!!! Diff         : ", m_mid.leash - mw.leash );
                println( "!!! lower_bound  : ", lower_bound );
                println( "!!! LB           : ", lb );
            end
            factor = factor * factor_scale;
            aprx_refinement = 1.0 + (aprx_refinement - 1.0) / approx_scale;
            f_debug  &&  println( "1 aprx_refinement: ", aprx_refinement,
                 "  scale :", approx_scale );

            pl, ql = Morphing_extract_vertex_radii( mw );
            continue;
        end

        f_debug  &&  println( "Extracting offsets..." );

        # Now we compute the distance, with offsets...
        #PS_offs = offsets_a;#Morphing_extract_offsets( m_a )[2]
        #QS_offs = offsets_b;#Morphing_extract_offsets( m_b )[1]

        f_debug  &&  println( "offsets_a max : ", maximum( offsets_a ) );
        f_debug  &&  println( "offsets_a max : ", maximum( offsets_b ) );

        f_debug  &&  println( "ve_r_mono( PS -> QS)" );

        map_PS = get_refinement_map( PSR, PS );
        map_QS = get_refinement_map( QSR, QS );

        PS_offs = offsets_copy_map( map_PS, offsets_a );
        QS_offs = offsets_copy_map( map_QS, offsets_b );
        @assert( length( PS_offs ) == cardin( PSR ) );
        @assert( length( QS_offs ) == cardin( QSR ) );

        # very important! We have to use the refined to monotonicity copies...
        m_final = frechet_ve_r_compute_ext( PSR, QSR, PS_offs, QS_offs,
                                            true );
        f_debug && println( "*** lower_bound.leash  : ", m_final.leash );
        f_debug && println( "*** mw     .leash      : ", mw.leash );
        #f_debug && println( "*** lower_bound       : ", lower_bound );

        if  ( eq( m_final.leash,  mw.leash, tolerance ) )
            f_debug && println( "Return from frechet_c_compute" );
            return  mw
        end

        factor = factor * factor_scale;
        aprx_refinement = 1.0 + ( (aprx_refinement - 1.0) /  approx_scale );
        f_debug  &&  println( "aprx_refinement: ", aprx_refinement,
             "  scale :", approx_scale );
        #lower_bound = lower_bound * 0.5; #max( lower_bound,
        #m_final.leash ); #lower_bound * 0.99;
        if f_debug
            println( "m_a.leash      : ", m_a.leash );
            println( "m_b.leash      : ", m_b.leash );
            println( "m_mid.leash    : ", m_mid.leash );
            println( "MW Leash length: ", mw.leash );
            println( "m_final.leash  : ", m_final.leash );
        end
        if  f_accept_approx  &&  ( 1.00001 * m_final.leash > mw.leash )
            f_debug && println( "MW Leash length: ", mw.leash );
            f_debug && println( "m_final.leash  : ", m_final.leash );
            f_debug && println( "Return (2) from frechet_c_compute" );
            return  mw
        end
    end  # While loop end
end

###########################################################################
###########################################################################
###########################################################################
# Greedy simplification


function  find_frechet_prefix_inner( P::Polygon{D,T},
    strt::Int64, i::Int64, j::Int64, w::T ) where {D,T}
    if  ( i >= j )  ||  ( (strt + 1) == j )
        return  j; # the width is zero, nothing to do.
    end
    r = frechet_width_approx( P, strt:j )
    if  ( r <= w )
        return  j;
    end
    if (i+1) == j
        return  i;
    end
    j = j - 1;
    mid = ( i + j ) >> 1;
    r_m = frechet_width_approx( P, strt:mid )
    if  ( r_m > w )
        return   find_frechet_prefix_inner( P, strt, i, mid - 1, w );
    end

    # No reason to waste time?
    if  (j - mid) <= ((j - strt) ÷ 10)
        #println( "FLOGI" );
        return  mid;
    end

    return  find_frechet_prefix_inner( P, strt, mid, j, w );
end

function  find_frechet_prefix( P::Polygon{D,T}, i::Int64,
    j::Int64, w::T )  where {D,T}
    r = frechet_width_approx( P, i:j )
    if  ( r <= w )
        return  j;
    end
    return  find_frechet_prefix_inner( P, i, i, j, w );
end

function  frechet_simplify_to_width( P::Polygon{D,T}, w::T ) where {D,T}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    len = cardin( P );
    if  ( len == 0 )
        return pout;
    end
    if  ( push_smart!( pout, P[1] ) )
        push!( pindices, 1 );
    end

    curr_ind = 1;
    next_ind::Int64 = 1;
    card = cardin( P );
    while  true
        next_ind = find_frechet_prefix( P, curr_ind, card, w )
        #println( "next_ind: ", next_ind );
        @assert( next_ind > curr_ind );
        push!( pindices, next_ind );
        push!( pout,  P[ next_ind ] );
        if  next_ind == card
            return  pout, pindices;
        end
        curr_ind = next_ind;
    end
end


function  frechet_simplify_to_cardin( P::Polygon{D,T}, sz::Int64 ) where {D,T}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    @assert( cardin( P ) >= 2 );
    if  sz <= 2
        push!( pindices, 1, cardin( P ) );
        return  spine( P ), pindices;
    end
    w = frechet_width_approx( P );
    ratio = 0.8;

    count::Int64 = 1;
    while  true
        count = count + 1;
        w = w * ratio;
        if  count > 5
            w = w / count;
        end
        Q, q_indices = Polygon_simplify_ext( P, w )
        m = frechet_c_mono_approx_subcurve( P, Q, q_indices )[ 1 ];
        R, R_indices = frechet_simplify_to_width( P, m.leash );
        if  ( ( cardin( R ) >= sz + 2 )  ||  ( 2*cardin( R ) >= cardin( P ) ) )
            println( " RRR :", cardin(R), " / ", sz );
            return  R, R_indices;
        end
        m = frechet_c_mono_approx_subcurve( P, R, R_indices )[ 1 ];
        w = m.leash;
    end
end


###########################################################################



function  exp_search_width_prefix( P::Polygon{N,T},
                                   start::Int64,
                                   w::Float64 ) where {N,T}
    len = cardin( P );
    hi::Int64 = min( start + 2, len );

    ( hi >= len )  &&  return  hi;
    while  hi < len
        r = frechet_width_approx( P, start:hi )
        ( r > w )  &&  return  min( hi + 5, len );
        hi = start + 2 * ( hi - start )
    end

    return  len;
end


function  frechet_simplify_w_exp( P::Polygon{D,T}, w::T ) where {D,T}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    card = cardin( P );
    if  ( card == 0 )
        return pout, pindices;
    end
    if  ( push_smart!( pout, P[1] ) )
        push!( pindices, 1 );
    end

    curr_ind = 1;
    next_ind::Int64 = 1;
    while  true
        hi = exp_search_width_prefix( P, curr_ind, w );
        next_ind = find_frechet_prefix( P, curr_ind, hi, w )
        @assert( next_ind > curr_ind );
        push!( pindices, next_ind );
        push!( pout,  P[ next_ind ] );
        if  next_ind == card
            return  pout, pindices;
        end
        curr_ind = next_ind;
    end
end


#
# End of file
##########################################################################

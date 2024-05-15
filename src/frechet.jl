# Originally contributed by S. Har-Peled
# under MIT License

#-------------------------------------------------------------
# Frechet
#
# A library to compute various versions of the Frechet distance.
#-------------------------------------------------------------


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

function  EID_i( id::Int64 )::Int64
    return  ( id >> 32 );
end

function  EID_j( id::Int64 )::Int64
    return   ( id & 0xfffffffc ) >> 2;
end

function   EID_i_is_vert( id::Int64 )::Bool
    return   ( (id & 0x2) != 0 );
end

function  EID_j_is_vert( id::Int64 )::Bool
    return  ( (id & 0x1) != 0 );
end

############################################################
# TreeVertex:
#    Encoding for a vertex in the VE Diagram how the search
#    procedure arrived ot it. I.e., encoding the reverse tree
#    that the algorithm computes.
#
# id: Vertex of the diagram
# val: The value associated with this vertex.
# id_prev: The vertex we arrvied from.
############################################################
@with_kw mutable struct  TreeVertex
    id::Int64
    val::Float64 = 0.0
    r::Float64 = 0.0

    id_prev::Int64
end

#DictVERType = Dict{Int64, TreeVertex};
DictVERType = Dict{Int64, Int64};
DictHandledType = Dict{Int64, Bool};
HeapVERType = BinaryMinHeap{TreeVertex};

function  TreeVertex( _id::Int64 )
    return TreeVertex( _id, 0.0, 0.0, EID( 0, false, 0, false ) );
end

function Base.isless( v::TreeVertex, u::TreeVertex )
    return  v.val < u.val;
end

###########################################################
# TODO: Seems unnecessary?
###########################################################
function Base.hash(v::TreeVertex, h::UInt)
    return  hash( v.id, hash( :TreeVertex, h) );
end


@with_kw mutable struct FRContext{N,T}
    P::Polygon{N,T};
    Q::Polygon{N,T};
    p_offs::Vector{Float64};
    q_offs::Vector{Float64};
    handled::DictHandledType;
    dict::DictVERType;
    heap::BinaryMinHeap{TreeVertex};
    f_offsets::Bool = false;
    n_p::Int64 = 0
    n_q::Int64 = 0
end

function  FRContext(P::Polygon{N,T}, Q::Polygon{N,T}) where {N,T}
    return FRContext( P, Q, Vector{Float64}(),  Vector{Float64}(),
                      Dict{Int64, Int64}(),
                      Dict{Int64, Float64}(),
                      BinaryMinHeap{TreeVertex}(),
                      false, cardin( P ), cardin( Q ) );
end

function  ve_event_value( c::FRContext{N,T}, id::Int64 ) where {N,T}

#    i::Int64, i_is_vert::Bool,
#                          j::Int64, j_is_vert::Bool ) where {N,T}
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
        d = dist_seg_nn_point( Q[ j ], Q[ j + 1 ], c.P[ i ] );
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
        d = dist_seg_nn_point( P[ i ], P[ i + 1 ], Q[ j ] );
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
    return  ev
end


function  f_r_new_event( _id::Int64, c::FRContext{N,T} ) where {N,T}
    ev = TreeVertex( _id );
    ev.r = ev.val = ve_event_value( c, ev.id );
    return  ev;
end

function   is_start_event( id::Int64 )
    return ( ( EID_i( id ) == 1 )
             && ( EID_j( id ) == 1 )
             && EID_i_is_vert( id )
             &&  EID_j_is_vert( id ) )
end

function  is_final_cell( id::Int64, n_p::Int64, n_q::Int64 )
    return   ( ( EID_i( id ) == ( n_p - 1 ) )
               &&  ( EID_j( id ) == ( n_q - 1 ) ) )
end

function f_r_create_event( R::Polygon{N,T}, i::Int64,
                        is_vert::Bool, qr::Point{N,T} ) where {N,T}
    if  ( is_vert )
        ev = EventPoint( R[ i ], i, PT_VERTEX, 0.0 );
        return  ev;
    end
    seg = Segment( R[ i ], R[ i + 1 ] );
    p = Segment_nn_point( seg, qr );

    t::Float64 = Segment_get_convex_coef( seg, p );

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

    return  EventPoint( deepcopy( p ), i, PT_ON_EDGE, t );
end


function  is_schedule_event( dict, id::Int64, n_p, n_q )::Bool
    if  haskey( dict, id )
        return  false
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

function  f_r_schedule_event( id::Int64, prev_id::Int64,
                              c::FRContext{N,T} ) where  {N,T}
    if  ! is_schedule_event( c.dict, id, c.n_p, c.n_q )
        return
    end
    ev = f_r_new_event( id, c);

    ev.id_prev = prev_id;
    c.dict[ id ] = prev_id;
    push!( c.heap, ev );
    return  ev
end




function  f_r_extract_solution( P::Polygon{N,T}, Q,
                                end_event_id::Int64, dict
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
    c::FRContext{N,T} = FRContext( P, Q )

    if  f_use_offsets
        c.p_offs = p_offs;
        c.q_offs = q_offs;
        c.f_offsets = true;
    end

    start_id = EID( 1, true, 1, true );

    end_id = EID( c.n_p, true, c.n_q, true );
    start_event = f_r_schedule_event( start_id, start_id, c );

    end_event::TreeVertex = start_event;
    iters = 0;
    heap = c.heap;
    while  ! isempty( heap )
        ev::TreeVertex = pop!( heap );
        if  ( haskey( c.handled,  ev.id ) )
            continue;
        end
        c.handled[ ev.id ] = true;
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

        value = ev.val;

        i = EID_i( ev.id );
        j = EID_j( ev.id );

        if  is_start_event( ev.id )
            f_r_schedule_event( EID( 1, false, 1, true ), ev.id, c );
            f_r_schedule_event( EID( 1, true, 1, false ), ev.id, c );
            continue;
        end

        # Is it the *final* event?
        if  ( i == c.n_p )  &&  ( j == c.n_q )
            end_event = ev;
            break;
        end

        # Is it on the boundary of the final cell?
        if  is_final_cell( ev.id, c.n_p, c.n_q )
            f_r_schedule_event( end_id, ev.id, c );
            continue;
        end
        f_r_schedule_event( EID( i+1, true, j, false ), ev.id, c );
        f_r_schedule_event( EID( i, false, j+1, true ), ev.id, c );
    end

    pes, qes = f_r_extract_solution( P, Q, end_event.id, c.dict );

    morph::Morphing{N,T} = Morphing_init( P, Q, pes, qes );
    morph.iters = iters;

    return  morph
end


function   frechet_ve_r_compute( P::Polygon{N,T}, Q::Polygon{N,T} ) where {N,T}
    return frechet_ve_r_compute_ext( P, Q, Vector{Float64}(),
                                     Vector{Float64}(), false );
end









#####################################################################
# 2-approximation to the Frechet distance between
#   P[first(rng)]-P[last(rng)] and he polygon
#   P[rng]
# Here, rng is a range i:j
#####################################################################
function  frechet_width_approx( P::Polygon{N,T},
                                rng::UnitRange{Int64} = 0:0
                                ) where {N,T}
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

    seg = Segment{N,T}( P[ first( rng ) ], P[ last( rng ) ] );

    t::Float64 = 0;
    curr::Point{N,T} = deepcopy( P[ first( rng ) ] ) ;
    leash::Float64 = 0;
    for  i  in  first(rng)+1:last(rng)-1
        q = Segment_nn_point( seg, P[ i ] );
        new_t = Segment_get_convex_coef( seg, q );
        if  ( new_t > t )
            t = new_t;
            curr = q;
        end
        leash = max( leash, Dist( curr, P[ i ] ) );
    end
    return  leash;
end

eachindexButLast(x) =  firstindex(x):lastindex(x)-1

#####################################################################
# frechet_offests

#    Given P, and subset of indices of P defining a subpolygon,
# computes a 2-approximation to the Frechet distance between every
# edge of the subpolygon they define, and corresponding portion of P.
#
#####################################################################
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
    lena = Polygon_length( polya );
    lenb = Polygon_length( polyb );

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
Q[ i ] =P[ p_indices[ i ].

"""
function    frechet_c_mono_approx_subcurve(
    P::Polygon{N,T},
    Q::Polygon{N,T},
    p_indices::Vector{Int64} ) where {N,T}

    @assert( cardin(Q) == length( p_indices ) );

    pes = Vector{EventPoint{N,T}}();
    qes = Vector{EventPoint{N,T}}();

    for  i  in  eachindexButLast( p_indices )
        pind = p_indices[ i ];
        pind_next = p_indices[ i + 1 ];
        push!( pes, EventPoint( P[ pind ], pind, PT_VERTEX, 0.0 ) );
        push!( qes, EventPoint( Q[ i ], i, PT_VERTEX, 0.0 ) );

        max_t::Float64 = 0;
        seg = Segment( Q[ i ], Q[ i + 1 ] );
        for  j in (pind + 1):(pind_next - 1)
            p = P[ j ];
            q = Segment_nn_point( seg, p );
            new_t = Segment_get_convex_coef( seg, q );
            if  ( new_t >= max_t )
                max_t = new_t;
            else # new_t < max_t
                q = Segment_get_on( seg, max_t )
            end
            push!( pes, EventPoint( P[ j ], j, PT_VERTEX, 0.0 ) );
            push!( qes, EventPoint( q, i, PT_ON_EDGE, max_t ) );
        end
    end

    push!( pes, EventPoint( last( P ), cardin( P ), PT_VERTEX, 0.0 ) );
    push!( qes, EventPoint( last( Q ), cardin(Q), PT_VERTEX, 0.0 ) );

    return  Morphing_init( P, Q, pes, qes );
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

    function  push_time( tm )
        if  ( tm == 0 )  ||  ( tm == 1 )
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

    for  tm in out_times
        if  ( tm == 0 ) || (tm == 1 )
            continue;
        end

        Polygon_push_smart( pout, Segment_get_on( seg, tm ) )
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
    extract_refined_polygon( P, pes, num_points_to_add )

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
                                   points_to_add::Int64 = 1 ) where {N,T}
    pout = Polygon{N,T}()
    opnts = Vector{EventPoint{N,T}}();

    len = length( s );
    push!( pout, s[1].p );

    i = 2;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
            Polygon_push_smart( pout, s[i].p );
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
        for  k in i:j
            push!( times, s[ k ].t );
        end

        if   is_monotone_inc( times )
            for  k in i:j
                Polygon_push_smart( pout, s[k].p );
            end
            i = j + 1
            continue
        end

        seg = Segment( poly[ loc ], poly[ loc + 1 ] );
        add_points_along_seg( pout, seg, times, points_to_add );

        i = j + 1
    end

    return  pout
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
        f_debug && println( "Before ve_r_compute( ", cardin(P), ", ",
                            cardin( Q ) );
        m = frechet_ve_r_compute( P, Q )
        f_debug  && println( "After... m.leash : ", m.leash );
        mm = Morphing_monotonize( m );
        f_debug && println( "mm.leash :", mm.leash );

        if  ( f_snapshots )
            push!( out, deepcopy( m ) );
        end

        fr_retract = m.leash;
        fr_r_mono = mm.leash;

        if ( f_debug )
            f_debug && println( "fr_retract : ", m.leash );
            f_debug && println( "fr_mono    : ", fr_r_mono );
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

    #println( "RETURNING!" );
    return  mm, f_exact, P, Q
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
function  frechet_ve_r_mono_compute( poly_a, poly_b )
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

Importantly, approx can be larger than 2, if you want a really
rough approximation.

"""
function  frechet_c_approx( poly_a::Polygon{N,T},
    poly_b::Polygon{N,T}, approx::Float64 ) where {N,T}

    @assert( ( cardin( poly_a ) > 1 )  &&  ( cardin( poly_b ) > 1 ) )
    f_debug::Bool = false;
    @assert( approx > 1.0 )

    f_do_one_round::Bool = true;
    d = frechet_dist_upper_bound( poly_a, poly_b )

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

            @assert( ( cardin( P ) > 1 )  &&  ( cardin( Q ) > 1 ) )

            f_debug &&  println( "before" );
            m = frechet_mono_via_refinement( P, Q,
                                             3.0/4.0 + approx / 4.0 )[1];
            d = m.leash;  #        frechet_ve_r_compute_dist( P, Q )
            f_debug && println( "after" );

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
        f_debug  &&  println( "f_c_dist( poly_a, P )" );
        m_p = frechet_c_mono_approx_subcurve( poly_a, P, p_indices );
#        println( "mono b" );

        #m_q = frechet_ve_r_mono_compute( Q, poly_b )
        f_debug  &&  println( "f_c_dist( poly_b, Q )" );
        m_q = frechet_c_mono_approx_subcurve( poly_b, Q, q_indices );
        Morphing_swap_sides!( m_q );

        err = max( m_p.leash, m_q.leash );
        #println( "Error: ", err );

        # m = m_p * mm * m_q
        f_debug && println( "Combine m_p, m" );
        mmu = Morphing_monotonize( Morphing_combine( m_p, m ) );
        f_debug && println( "Combine mmu,m_q" );
        mm_out = Morphing_combine( mmu, m_q );

        f_debug && println( "Verification..." );
        Morphing_verify_valid( m_p );
        Morphing_verify_valid( m );
        Morphing_verify_valid( m_q );
        Morphing_verify_valid( mmu );

        ratio = mm_out.leash/(d - 2.0*err);
        f_debug &&  println( "Ratio :", ratio );
        if  ratio <=  approx
            if f_debug
                println( "-----------------------------------" );
                println( "Leash out: ", mm_out.leash );
                println( "Leash lower bound: ", d - 2r );
                println( "approx ratio          : ",ratio  );
                println( "Required approximation: ", 1.0 + approx );
                println( "||P||: ", Polygon_length( poly_a ),
                    "  ||Q||: ",
                    Polygon_length( poly_b ) );
            end
            mm_out.ratio = ratio;
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

    PS, p_indices = Polygon_simplify_radii( P, r );
    if   ( cardin( PS ) <= ( cardin( P ) / 2 ) )
        return  PS, p_indices, false
    end

    pindices = Vector{Int64}();
    for  i in cardin( P )
        push!( pindices, i )
        return  P, p_indices, true
    end
end

mutable struct FrechetCExtraInfo
    PSR::Polygon2F;
    QSR::Polygon2F;
    PSR_offs::Vector{Float64};
    QSR_offs::Vector{Float64};
    f_init::Bool
end;

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
function  frechet_c_compute( poly_a::Polygon{N,T},
                             poly_b::Polygon{N,T},
                             f_accept_approx::Bool = true
                             )  where {N,T}
    f_debug::Bool = false;
    aprx_refinement::Float64 = 1.1;

    f_debug  &&  println( "\n\n\n\n\n\n" );
    f_debug  &&  println( "#", cardin( poly_a ) )
    f_debug  &&  println( "#", cardin( poly_b ) )

    mf = frechet_c_approx( poly_a, poly_b, 2.0 );
    ratio_2 = mf.ratio;

    len_a = Polygon_length( poly_a );
    len_b = Polygon_length( poly_b );
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
        ratio = ratio_2;
    else
        f_debug  &&  println( "freceht_c_approx( ", approx, ") " );
        m = frechet_c_approx( poly_a, poly_b, approx );
        ratio = m.ratio
        f_debug  &&  println( "freceht_c_approx( ", approx, ")...done" );
    end
        #        frechet_ve_r_mono_approx( poly_a, poly_b, eps );

    if  ( f_debug )
        println( "Approximation computed..." );
        println( "m.leash : ", m.leash );
        println( "" );
    end

    pl, ql = Morphing_extract_vertex_radii( m );

    lower_bound = m.leash / ratio;
    if  f_debug
        println( "Before while loop..." );
        println( "LOWER BOUND: ", lower_bound );
        println( "upper BOUND: ", m.leash );
    end
    f_debug  &&  println( "\n\n\n\n\n\n" );
    factor::Float64 = 64.0
    while  true
        f_debug  &&  println( "-------------------------------------------" );
        f_debug  &&  println( "factor: ", factor, "                      " );
        pz = ( ( lower_bound * ones( length( pl ) ) ) - pl ) / factor
        qz = ( ( lower_bound * ones( length( ql ) ) ) - ql ) / factor

#        f_deprintln( pz );

        p_count = count_below_zero( pz );
        q_count = count_below_zero( qz );


        f_PS_exact::Bool = false;
        f_QS_exact::Bool = false;

        PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext( poly_a, pz );
        QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext( poly_b, qz );

        if  f_debug
            println( "PS.len    : ", cardin( PS ); );
            println( "QS.len    : ", cardin( QS ); );
            println( "p_count <0: ", p_count, " / ", length( pz ) );
            println( "q_count <0: ", q_count, " / ", length( qz ) );
            println( "|PS| = ", cardin( PS ) );
            println( "|P| = ", cardin( poly_a ) );
            println( "|QS| = ", cardin( QS ) );
            println( "|Q| = ", cardin( poly_b ) );
            println( "Computing radii simplified Frechet distance..." );
        end
        #    m_mid = frechet_ve_r_mono_compute( PS, QS  );
        f_debug && println( "\nApprox refinement : ", aprx_refinement );
        m_mid, f_exact, PSR, QSR = frechet_mono_via_refinement( PS, QS,
                                                              aprx_refinement );

        f_debug  &&  println( "frechet mono via refinment computed" );
        f_debug  &&  println( "PSR.len: ", cardin( PSR ) );
        f_debug  &&  println("QSR.len: ", cardin( QSR ) );

        f_debug  &&  println( "ve_r_mono( poly_a -> PSR)" );
        m_a = frechet_ve_r_mono_compute( poly_a, PSR );
        mmu = Morphing_combine( m_a, m_mid );
        f_debug  &&  println( "ve_r_mono( WSR -> poly_b )" );
        m_b = frechet_ve_r_mono_compute( QSR, poly_b );
        mw = Morphing_combine( mmu, m_b );

        # is there is hope we got the optimal solution?
        if  ( ! floating_equal( m_mid.leash, mw.leash ) )
            if  f_debug
                println( "m_mid.leash  : ", m_mid.leash );
                println( "mw.leash     : ", mw.leash );
                println( "Diff         : ", m_mid.leash - mw.leash );
            end
            factor = factor * 2.0;
            aprx_refinement = 1.0 + (aprx_refinement - 1.0) / 4.0;
            continue;
        end

        # Now we compute the distance, with offsets...
        PSR_offs = Morphing_extract_offsets( m_a )[2]
        QSR_offs = Morphing_extract_offsets( m_b )[1]


        f_debug  &&  println( "ve_r_mono( PSR -> QSR)" );
        m_final = frechet_ve_r_compute_ext( PSR, QSR, PSR_offs, QSR_offs,
                                            true );
        if  ( floating_equal( m_final.leash,  mw.leash ) )
            f_debug && println( "Return from frechet_c_compute" );
            return  mw
        end

        factor = factor * 2.0;
        aprx_refinement = (1.0 + aprx_refinement) / 2.0;
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


#
# End of file
##########################################################################

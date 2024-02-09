##############################################################
# Frechet
#
# A library to compute various versions of the Frechet distance.
###############################################################

using Parameters
using DataStructures
using Printf

#push!(LOAD_PATH, pwd())
#using cg

#########################################################################3
############################################################################
# Retractable Frechet distance
############################################################################

@enum FPointType begin
    PT_VERTEX = 1
    PT_ON_EDGE = 2
end

############################################################
# A vertex edge event descriptor.
#
# p: Location of the point being matched. Not necessarily a vetex of
#    the polygon.
# i: Vertex number in polygon/or edge number where p lies.
# type:
# t: Convex combination coefficient if p is on the edge. t=0 means its
#    the ith vertex, t=1 means it is on the i+1 vertex.
############################################################
#@with_kw
mutable struct EventPoint{N,T}
    p::Point{N,T}
    i::Int64;
    type::FPointType;
    t::Float64;
end

# leash = maximum length of an edge in the encoded matching.
#         (That is, the Frechet distance.)
"""
    Morphing

Encoding of a moprhing (i.e., matching) between two polygonal cuves.

"""
@with_kw mutable struct  Morphing{N,T}
    P::Polygon{N,T}
    Q::Polygon{N,T}
    pes::Vector{EventPoint{N,T}}; # P event sequence
    qes::Vector{EventPoint{N,T}}; # Q event sequence
    leash::Float64
    leash_offsets::Float64
    iters::Int64
    ratio::Float64
end

function   Morphing_init( P::Polygon{N,T}, Q::Polygon{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}} ) where {N,T}

    @assert( length( pes ) == length( qes ) );
    m = Morphing( P, Q, pes, qes, 0.0, 0.0, 0, 0.0 );
    Morphing_recompute_leash( m );

    return  m;
end


####################################################
# Compute for every vertex the maximum leash length used with this
# vertex.
####################################################

function   extract_vertex_leash( P::Polygon{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}}
) where {N,T}
    P_len = cardin( P );
    len = length( pes );
    vl = zeros( Float64, P_len )
    for  i in  1:len
        ev = pes[ i ];
        evx = qes[ i ];
        if  ( ev.type == PT_VERTEX )
            vl[ ev.i ] = max( vl[ ev.i ], Dist( ev.p, evx.p ) );
        end
    end

    return  vl;
end

function   Morphing_extract_vertex_radii( m::Morphing{N,T} ) where {N,T}
    pl = extract_vertex_leash( m.P, m.pes, m.qes );
    ql = extract_vertex_leash( m.Q, m.qes, m.pes );
    return  pl, ql
end


function   Morphing_swap_sides!( m::Morphing{N,T} ) where {N,T}
    m.P, m.Q = m.Q, m.P
#    swap!( m.P, m.Q );
    #swap!( m.pes, m.qes );
    m.pes, m.qes = m.qes, m.pes;
end

function   Morphing_empty( P::Polygon{N,T}, Q::Polygon{N,T} )  where {N,T}
    pes = Vector{EventPoint{N,T}}();
    qes = Vector{EventPoint{N,T}}()
    r::Float64 = -1;
    return  Morphing( P, Q, pes, qes, r, 0.0, 0, 0.0 );
end


###################################################################
# Check if the matching is monotone. That is the endpoints of the
# matching on each curve are sorted, in the same order as the
# matching.
##################################################################
function  events_seq_is_monotone( s::Vector{EventPoint{N,T}}
) where {N,T}

    len = length( s );
    i::Int64 = 1;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
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

        # i:j is the sequence of edge-vertex events
        for  k  in i:(j-1)
            if  s[ k + 1 ].t < s[ k ].t
	        return  false;
            end
        end
        i = j + 1;
    end

    return  true;
end


##############################################################################
# Make the matching (which might be potentially not monotone) into a monotone
# one.
##############################################################################
function  events_seq_make_monotone( P::Polygon{N,T},
                                    s::Vector{EventPoint{N,T}} ) where {N,T}
    ns = Vector{EventPoint{N,T}}();

    len = length( s );
    i::Int64 = 1;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
            push!( ns, deepcopy( ep ) )
            i = i + 1;
            continue;
        end

        j = i
        loc = s[ i ].i;
        t = s[ i ].t
        while  ( ( j < len )
                 &&  ( s[ j + 1 ].type == PT_ON_EDGE )
                 &&  ( s[ j + 1 ].i == loc ) )
            nep = deepcopy( s[ j ] );
            t = max( t, nep.t );
            if  ( t > nep.t )
                nep.t = t;
                ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
                nep.p = convex_comb( P[ nep.i ], P[ nep.i + 1 ], t );
            end
            push!( ns, nep );
            j = j + 1
        end

        nep = deepcopy( s[ j ] );
        if  ( t > nep.t )
            nep.t = t;
            nep.p = convex_comb( P[ nep.i ], P[ nep.i + 1 ], t );
            ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
        end
        push!( ns, nep );

        i = j + 1;
    end

    @assert( length( ns ) == length( s ) );
    return  ns;
end

function  Morphing_recompute_leash( m::Morphing{N,T} ) where  {N,T}
    @assert( length( m.pes ) == length( m.qes ) );
    r = 0;
    for  i  in eachindex( m.pes )
        ell = Dist( m.pes[ i ].p, m.qes[ i ].p )
        if   ell > r
            r = ell;
        end
    end
    if  ( m.leash != r )   &&   ( m.leash != 0 )
        println( "old leash:" , m.leash );
        println( "new leash:" , r );
        exit( -1 );
    end
    m.leash = r;
#    println(  " m.leash: ", m.leash );
end

######################################################################
# Turng the morphing matching into a "true" matching, by creating
# two polygons that their edges are directly matched.
# The output polygons P and Q will definitely have reapeated points.
######################################################################
function  Morphing_as_polygons( m::Morphing{N,T} ) where  {N,T}
    P = Polygon{N,T}();
    Q = Polygon{N,T}();

    for  i  in eachindex( m.pes )
        push!( P, m.pes[ i ].p );
        push!( Q, m.qes[ i ].p );
    end

    return  P, Q
end


function   Morphing_is_monotone( m::Morphing{N,T} ) where {N,T}
    return ( events_seq_is_monotone( m.pes )
             && events_seq_is_monotone( m.qes ) )
end

##########################################################3
# Turn the morphing into monotone morphing...
##########################################################3
"""
    Morphing_monotonize

Turns a morphing into a monotone morphing, by simply not going back,
staying in place if necessary.

"""

function  Morphing_monotonize( m::Morphing{N,T} ) where {N,T}
    if  Morphing_is_monotone( m )
        return  deepcopy( m );
    end
    P = m.P;
    Q = m.Q;

    pes_new = events_seq_make_monotone( P, m.pes );
    qes_new = events_seq_make_monotone( Q, m.qes );
    return  Morphing_init( P, Q, pes_new, qes_new );
end



##########################################################################
##########################################################################
# Computing the VE Frechet distance...
##########################################################################
##########################################################################

# EID = vertices of the VE Frechet diagram.
@with_kw mutable struct  EIDOLD
    i::Int64 = 0;
    i_is_vert::Bool = false;
    j::Int64 = 0;
    j_is_vert::Bool = false;
end

function  EID(
    i::Int64 = 0,
    i_is_vert::Bool = false,
    j::Int64 = 0,
    j_is_vert::Bool = false
)
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
    dict::DictVERType;
    heap::BinaryMinHeap{TreeVertex};
    f_offsets::Bool = false;
    n_p::Int64 = 0
    n_q::Int64 = 0
end

function  FRContext(P::Polygon{N,T}, Q::Polygon{N,T}) where {N,T}
    return FRContext( P, Q, Vector{Float64}(),  Vector{Float64}(),
                      Dict{Int64, Int64}(),
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
            return   Dist( P[ i ], Q[ j ] );
        end
        d = dist_seg_nn_point( Q[ j ], Q[ j + 1 ], c.P[ i ] );
        if  c.f_offsets
#            if  ( c.q_offs[ j ] > 0 )
#                println( "Larger than 0!" );
#            end
            return  d - c.q_offs[ j ];
        else
            return  d;
        end
    end

    if  EID_j_is_vert( id )
        d = dist_seg_nn_point( P[ i ], P[ i + 1 ], Q[ j ] );
        if  c.f_offsets
            return  d - c.p_offs[ i ];
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
    @assert( 0.0 <= t  &&  t <= 1.0 );
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

function  print_int_w_commas( n::Int64 )
    if  ( abs( n ) < 1000 )
        print( n )
        return
    end
    print_int_w_commas( floor( Int64, n/1000 ) );
    @printf( ",%03d",  n % 1000 )
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
                                     f_use_offsets::Bool = false
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
        iters = iters + 1;

        if  f_debug  &&  ( (iters % 1000000) == 0 )
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



############################################################################
############################################################################
# Discrete frechet
#--------------------------------------------------------------------------


function   d_frechet_extract_solution( P::Polygon{N,T}, Q,
    dp_dec_i, n_p, n_q )::Morphing{N,T} where {N,T}
    i = n_p;
    j = n_q;

    println( "Extract solution?" );
    up::Int64 = n_p + n_q + 1;
    peout = Vector{EventPoint{N,T}}( undef, up );
    qeout = Vector{EventPoint{N,T}}( undef, up );

    pos::Int64 = 0;
    while  ( ( i != 1 )  ||  ( j != 1 )  )
        pos = pos + 1;
        peout[ pos ] = EventPoint( P[ i ], i, PT_VERTEX, 0.0 );
        qeout[ pos ] = EventPoint( Q[ j ], j, PT_VERTEX, 0.0 );

        if  dp_dec_i[ i, j ]
            i = i - 1;
        else
            j = j - 1;
        end
    end
    pos = pos + 1;
    peout[ pos ] = EventPoint( P[ 1 ], i, PT_VERTEX, 0.0 );
    qeout[ pos ] = EventPoint( Q[ 1 ], j, PT_VERTEX, 0.0 );

    @assert( pos <= up );
    resize!( peout, pos );
    resize!( qeout, pos );

    reverse!( peout );
    reverse!( qeout );

    return  Morphing_init( P, Q, peout, qeout );
end



function   frechet_d_compute_inner( P::Polygon{N,T}, Q::Polygon{N,T},
                                    dp::Array{Float64, 2}, dp_dec_i
                                        ) where {N,T}
    d::Float64 = 0;
    jp::Int64 = 0;
    ip::Int64 = 0;

    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    dp[ 1, 1 ] = Dist( P[ 1 ], Q[ 1 ] );
    for  i::Int64 in 1:n_p
        ip = max( i - 1, 1 );
        for  j in 1:n_q
            d = DistSq( P[ i ], Q[ j ] );
            jp = max( j - 1, 1 );

            f_dec_i = false;
            if  ( j == 1 )
                f_dec_i = true;
            else
                if  ( i == 1 )
                    f_dec_i = false;
                else
                    if  ( dp[ i - 1, j ] < dp[ i, j - 1 ] )
                        f_dec_i = true;
                    else
                        f_dec_i = false;
                    end
                end
            end

            if  ( f_dec_i )
                dp[ i, j ] = max( d, dp[ ip, j ] );
                dp_dec_i[ i, j ] = true;
            else
                dp[ i, j ] = max( d, dp[ i, jp ] );
                dp_dec_i[ i, j ] = false;  # redundant
            end
        end
    end

#    return  1.0
end

function   frechet_d_compute( P::Polygon{N,T},
                              Q::Polygon{N,T} ) where {N,T}
    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    dp::Array{Float64, 2} = Array{Float64, 2}(undef, n_p,n_q);
    dp_dec_i = falses( n_p, n_q );

    frechet_d_compute_inner( P, Q, dp, dp_dec_i );

    return  d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
end


function   frechet_d_compute_dist( P::Polygon{N,T},
                                   Q::Polygon{N,T} ) where {N,T}
    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    dp::Array{Float64, 2} = Array{Float64, 2}(undef, n_p,n_q);
    dp_dec_i = falses( n_p, n_q );

    frechet_d_compute_inner( P, Q, dp, dp_dec_i );

    return  dp[ n_p, n_q ];
end



##########################################################################
# Compute discrete frechet distance that is locally optimal
# Frechet. Essentially discrete frechet + Prim/Dijkstra algorithm For
# the discrete case!
##########################################################################
# _lopt_frechet
function   frechet_d_r_compute( P::Polygon{N,T}, Q::Polygon{N,T}
                                            ) where {N,T}
    n_p = cardin( P );
    n_q = cardin( Q );

    dp = zeros( n_p, n_q );
    dp_dec_i = falses( n_p, n_q );
    dp_def = falses( n_p, n_q );
    in_heap = falses( n_p, n_q );

    dp[ 1, 1 ] = Dist( P[ 1 ] , Q[ 1 ] );
    dp_def[ 1, 1 ] = true;

    x = [( Float64(dp[ 1, 1 ]), (1, 1) ) ];

    heap = BinaryMinHeap( x );
    in_heap[ 1, 1 ] = true;

    iters = 0;
    while  ! isempty( heap )
        ele = pop!( heap );
        iters = iters + 1;

        coords = ele[2];
        i = coords[ 1 ];
        j = coords[ 2 ];
        value = ele[ 1 ];

        d = Dist( P[ i ], Q[ j ] );

        ip = max( i - 1, 1 );
        jp = max( j - 1, 1 );

        f_dec_i::Bool = false;
        if      ( j == 1 )
            f_dec_i = true;
        elseif  ( i == 1 )
            f_dec_i = false;
        elseif  ! dp_def[ i - 1, j ]
            f_dec_i = false;
        elseif  ! dp_def[ i, j - 1 ]
            f_dec_i = true;
        else
            if  ( dp[ i - 1, j ] < dp[ i, j - 1 ] )
                f_dec_i = true;
            else
                f_dec_i = false;
            end
        end

        if  ( f_dec_i )
            dp[ i, j ] = max( d, dp[ ip, j ] );
            dp_def[ i, j ] = true;
            dp_dec_i[ i, j ] = true;
        else
            dp[ i, j ] = max( d, dp[ i, jp ] );
            dp_def[ i, j ] = true;
            dp_dec_i[ i, j ] = false;  # redundant
        end

        # Now we need to schedule the next two adjacent entries..
        inext = i + 1;
        jnext = j + 1;
        if  ( inext <= n_p )
            if  ( ! in_heap[ inext, j ] )
                push!( heap, ( Dist( P[ inext ], Q[ j ] ),
                               (inext, j) ) );
                in_heap[ inext, j ] = true;
            end
        end
        if  ( jnext <= n_q )
            if  ( ! in_heap[ i, jnext ] )
                push!( heap, ( Dist( P[ i ], Q[ jnext ] ),
                               (i, jnext) ) );
                in_heap[ i, jnext ] = true;
            end
        end
        if  ( ( i == n_p )  &&  ( j == n_q ) )
            break;
        end
    end

    println( "Iters: ", iters, " : ", n_p * n_q );
    return   d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
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

function  check_times( V::Vector{EventPoint{N,T}} ) where {N,T}
    for  ev in V
        if  ( ev.t < 0.0 )  ||  (ev.t > 1.0 )
            println( "Event time is wrong? ", ev.t );
            exit( -1 );
        end
    end
end


# Check the times stemps of the events are valid.
function  Morphing_verify_valid( m::Morphing{N,T} ) where {N,T}
    check_times( m.pes );
    check_times( m.qes );
end


"""
    extract_param_inner

Takes the sequence of events (i.e., sequence of points along the
curve), and outputs a corresponding sequence of real numbers, where
the ith nubmer if the distance of the ith point in sequence from the
beginning of the curve. Here, distance means the total distance of the
subcurve from the start point, to the current point.
"""
function  extract_param_inner( P::Polygon{N,T},
                               V::Vector{EventPoint{N,T}} ) where {N,T}
    lens::Vector{Float64} = Polygon_prefix_lengths( P )
    out = Vector{Float64}();
    check_times( V );

    n = cardin( P );
    for  i in eachindex( V )
        ev::EventPoint{N,T} = V[ i ];
        if  ev.type == PT_VERTEX
            push!( out, lens[ ev.i ] );
            continue;
        end
        #println( "ev.t: ", ev.t, "  i: ", i,"  len(V): ", length( V )  );
        curr = lens[ ev.i ];
        next = lens[ min( n, ev.i + 1 ) ]
        push!( out, curr + ev.t * (next - curr) )
    end

    return  out
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
    Morphing_extract_prm

    A parameterization is a polygonal curve that starts at (0,0) and
    end at (m,n).  The polygonal curve ither have positive slope edge,
    or vertical or horizontla edges. It can be thought of as a
    piecewise linear function from [0,m] to [0,n]. Here m and n are the
    lengths of the two given polygons of P and Q, respectively.

"""
function  Morphing_extract_prm( m::Morphing{N,T} )::Polygon2F where {N,T}
    P::Polygon{N,T} = m.P;
    Q::Polygon{N,T} = m.Q;
    peout::Vector{EventPoint{N,T}} = m.pes;
    qeout::Vector{EventPoint{N,T}} = m.qes;

    check_times( m.pes );

    check_times( m.qes );

    # pps and qps are two
    pps::Vector{Float64} = extract_param_inner( P, peout );
    qps::Vector{Float64} = extract_param_inner( Q, qeout );

    @assert( length( pps ) == length( qps ) )

    out = Polygon{2,Float64}();
    for  i in eachindex(pps)
        push!( out, point( pps[ i ], qps[ i ] ) )
    end
    return  out;
end


function  eval_pl_func_on_dim( p::Point{N,T}, q::Point{N,T}, val::Float64,
                               d::Int64 )  where {N,T}
    @assert( p[ d ] <= q[ d ] );
    t = (val - p[d]) / (q[d] - p[d]);

    return  p * (1.0 - t) + q * t;
end

function  eval_pl_func( p::Point2F, q::Point2F, val::Float64 )
    o = eval_pl_func_on_dim( p, q, val, 1 );
    return  o[2];
end
function  eval_inv_pl_func( p::Point2F, q::Point2F, val::Float64 )
    o = eval_pl_func_on_dim( p, q, val, 2 );
    return  o[1];
end

function  check_monotone_top( out::Polygon2F )
    l = cardin( out )
    if  l < 2
        return;
    end
    p = out[ l - 1 ];
    q = out[ l ];
    factor::Float64 = 1.0001
    if  ( p[ 1 ] > factor*q[1] ) || ( p[2] > factor*q[2] )
        println( p[ 1 ], "  >  ", q[1] );
        println( p[ 2 ], "  >  ", q[2] );
        println( "Error (not monotone top): " );
        for i in 1:l
            println( i, " : ", out[ i ] );
        end
        #println( out );
        println( "Error (not monotone top): " );
        exit( -1 );
    end
end

function  floating_equal( a::Float64, b::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (0.0000001* (abs( a)  + abs(b) ))
end


function  floating_equal( a, b )::Bool
    return  a == b;
end

##############################################################3
# parameterization_combine
# Given two parameterization f, g, compute f(g( . ))
###############################################################
function   parameterization_combine( f::Polygon2F, g::Polygon2F )
    if  ( ! floating_equal( last( g )[2], last( f )[ 1 ] ) )
        println( last( g )[2], " != ",  last( f )[ 1 ] );
        @assert( floating_equal( last( g )[2], last( f )[ 1 ] ) )
    end
    idf::Int64 = 1;
    idg::Int64 = 1;
    out::Polygon2F = Polygon2F();

    l_f = cardin( f )
    l_g = cardin( g )

    # x = dim(g,1)
    # y = dim(g,2) = dim(f, 1)
    # z = dim(f,2)
    while (true)
        check_monotone_top( out );

        # The two points under consideration, on the y axis
        yf = f[ idf ][1];
        yg = g[ idg ][2];
        f_equal_yf_yg::Bool = floating_equal( yf, yg )
        if  ( f_equal_yf_yg )  &&  ( idf == l_f )  &&  ( idg == l_g )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            break;
        end
        if  ( f_equal_yf_yg   &&  ( idf < l_f )
              &&  ( f[ idf + 1 ][1] == yf ) )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idf = idf + 1
            continue;
        end
        if  ( f_equal_yf_yg   &&  ( idg < l_g )
              &&  ( floating_equal( g[ idg + 1 ][2], yg ) ) )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idg = idg + 1
            continue;
        end
        if  f_equal_yf_yg
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idf = min( idf + 1, l_f );
            idg = min( idg + 1, l_g );
            continue;
        end
        if  ( yf < yg )
            # compute g( yf )...
            xg = eval_inv_pl_func( g[idg - 1], g[ idg ], yf )
    #        println( "xg: ", xg );

            # A bit of a hack... Because of floating point errors,
            # things can go a bit backward, which we really should not
            # allow.
            if  ( xg < g[ idg - 1 ][ 1 ] )
                xg = g[ idg - 1 ][ 1 ]
            end
            push!( out, point( xg, f[ idf ][ 2 ] ) )
            idf = min( idf + 1, l_f );
            continue;
        end
        if  ( yf > yg )
            zf = eval_pl_func( f[idf - 1], f[ idf ], yg )
            push!( out, point( g[ idg ][ 1 ], zf ) )
            idg = min( idg + 1, l_g );
            continue;
        end
        @assert( false );
    end

    # We monotononize the output because of floating point error...
    # Should really just get rid of tiny flatting point jitters.
    curr::Point2F = deepcopy( out[ 1 ] );
    for  i in 2:cardin( out )
        for  j in 1:2
            curr[ j ] = max( curr[ j ], out[ i ][ j ] )
            out[ i ][ j ] = curr[ j ];
        end
    end

    return  out;
end


function  get_point(
    P::Polygon{N,T}, lens::Vector{Float64},
    i::Int64, pos::Float64 ) where {N,T}

    if  ( i == length( lens ) )
        return  last( P ),1;
    end
    edge_len = lens[ i + 1 ] - lens[ i ];
    t = ( pos - lens[ i ] ) / edge_len;

    if  ( 0 > t > -0.00000001 )
        t = 0;
    end
    if  ( 1.0 < t < 1.000000001 )
        t = 1;
    end

#    println( "inside... 3" );
    if  ( t == 0 )
        return  P[ i ], .00;
    end
    if  ( t == 1 )
        return  P[ i + 1 ], 1.0;
    end
    p =  P[ i ] * (1.0 - t) + P[i + 1] * t;
#    println( "inside... 4" );
    return p, t;
end

function  event_sequences_extract( prm::Polygon2F, P::Polygon{N,T},
                                   Q::Polygon{N,T} ) where {N,T}

    i_p::Int64 = 1;
    i_q::Int64 = 1;
    pes = Vector{EventPoint{N,T}}(); # P event sequence
    qes = Vector{EventPoint{N,T}}(); # Q event sequence

    len_q = cardin( Q );
    len_p = cardin( P );

    lp::Vector{Float64} = Polygon_prefix_lengths( P )
    lq::Vector{Float64} = Polygon_prefix_lengths( Q )

    for  i in 1:( cardin( prm ) - 1 )
#        println( "before?? A" );
        curr::Point2F = prm[ i ];
        p_loc = curr[ 1 ];

        while  ( i_p < length( lp ) )  &&  ( p_loc >= lp[ i_p + 1 ] )
            i_p = i_p + 1;
        end
        if  ( ! ( lp[ i_p ] <= p_loc ) )
            println( "-----------------------" );
            println( "lp[ i_p - 1]  : ", lp[ i_p - 1 ] );
            println( "lp[ i_p ]     : ", lp[ i_p ] );
            println( "lp[ i_p + 1 ] : ", lp[ i_p + 1] );
            println( "p_loc         : ", p_loc );
            println( " i_p : ", i_p );
            println( " len : ", length( lp ) );
            println( prm[ i - 1 ] );
            println( prm[ i     ] );
            println( prm[ i + 1 ] );
#            println( prm );

        end
        #println( lp[ i_p ], "<=", p_loc );
        @assert( lp[ i_p ] <= p_loc  );
        while  ( i_q < length( lq ) )  &&  ( curr[ 2 ] >= lq[ i_q + 1 ] )
            i_q = i_q + 1;
        end

        #        println( "get_point P..." );
        pcurr,t_p = get_point( P, lp, i_p, p_loc );
        @assert( 0 <= t_p  &&  t_p <= 1.0 )
        #        println( "get_point Q ..." );
        get_point( Q, lq, i_q, curr[ 2 ] );
#        println( "get_point Q!!!!!!!!!! ..." );
        qcurr,t_q = get_point( Q, lq, i_q, curr[ 2 ] );

        if  ( ! ( 0.0 <= t_q <= 1.0 ) )
            println( "t_q: ", t_q );
            @assert( 0.0 <= t_q <= 1.0 );
        end
        #        println( "here? " );

        if  ( t_p == 0.0 ) ||  ( ( t_p == 1.0 )  &&  ( i_p == len_p ) )
            push!( pes, EventPoint( pcurr, i_p, PT_VERTEX, 0.0 ) );
        else
            push!( pes, EventPoint( pcurr, i_p, PT_ON_EDGE, t_p ) );
        end

        if  ( t_q == 0.0 ) ||  ( ( t_q == 1 )  &&  ( i_q == len_q ) )
            push!( qes, EventPoint( qcurr, i_q, PT_VERTEX, 0.0 ) );
        else
            push!( qes, EventPoint( qcurr, i_q, PT_ON_EDGE, t_q ) );
        end
        #       println( "here? B" );
    end

    push!( pes, EventPoint( last(P), cardin( P ), PT_VERTEX, 0.0 ) );
    push!( qes, EventPoint( last(Q), cardin( Q ), PT_VERTEX, 0.0 ) );

    #=
    for  i in 1:length( qes )
        if ( ( pes[ i ].type == PT_ON_EDGE ) &&
            ( qes[ i ].type != PT_ON_EDGE ) )
            println( "double edge event?" );
            println( "i:", i );
            println( "qes[ i ].t:", qes[ i ].t );
            println( "pes[ i ].t:", pes[ i ].t );
            println( length( qes ) );
        end
        @assert( ( pes[ i ].type != PT_ON_EDGE )
                 ||  ( qes[ i ].type != PT_ON_EDGE ) );
    end
    =#
    return  pes, qes
end

########################################################################
"""
    Morphing_from_prm( prm, P, Q )

    Construct a morphing from a parameterization of the two polygons P
    and Q.

"""
function  Morphing_from_prm( prm::Polygon2F,
                             P::Polygon{N,T},
                             Q::Polygon{N,T} ) where {N,T}
    pes, qes = event_sequences_extract( prm, P, Q );
    return  Morphing_init( P, Q, pes, qes );
end


"""
    add_points_along_seg(

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



function   compute_offsets( P::Polygon{N,T},
                            pes::Vector{EventPoint{N,T}},
                            qes::Vector{EventPoint{N,T}},
                            )  where  {N,T}
    offs = zeros( Float64, cardin( P ) );

    for  i  in eachindex( pes )
        ev = pes[ i ];
        evx = qes[ i ];
        r = Dist( ev.p, evx.p );
        offs[ ev.i ] = max( offs[ ev.i ], r )
    end

    return  offs
end

"""
    Morphing_extract_offsets
"""
function  Morphing_extract_offsets( m::Morphing )
    p_offs = Vector{Float64}();
    q_offs = Vector{Float64}();


    p_offs = compute_offsets( m.P, m.pes, m.qes );
    q_offs = compute_offsets( m.Q, m.qes, m.pes );
    ; # P event sequence

    return  p_offs, q_offs
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
"""
function  frechet_mono_via_refinement( Pa::Polygon{N,T}, Qa::Polygon{N,T},
                                       approx::Float64 = 1.00001 )  where {N,T}
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



###########################################################################
# Compute the morphing src(v) -> trg(u): That is u(v(t))
###########################################################################

"""
    Morphing_combine

Gets two morphings u, v (i.e., two parameterizations) and combine them
into a single morphing u(v(.)).

For example, if u: γ → δ  and  v: δ → ξ, then the returned morphing is
u(v(⋅)): γ → ξ.

"""
function  Morphing_combine( u::Morphing{N,T}, v::Morphing{N,T} ) where {N,T}
    u_prm = Morphing_extract_prm( u );
    v_prm = Morphing_extract_prm( v );

    #println( "b5678" );
    prm = parameterization_combine( v_prm, u_prm );

    pes, qes = event_sequences_extract( prm, u.P, v.Q );

    m = Morphing_init( u.P, v.Q, pes, qes );

    return  m;
end

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
                             f_accept_approx::Bool = true )  where {N,T}
    f_debug::Bool = false;
    aprx_refinement::Float64 = 1.001;

    f_debug && println( "#", cardin( poly_a ) )
    f_debug && println( "#", cardin( poly_b ) )

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
    factor::Float64 = 4.0
    while  true
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
            println( "p_count <0: ", p_count, " / ", length( pz ) );
            println( "q_count <0: ", q_count, " / ", length( qz ) );
            println( "|PS| = ", cardin( PS ) );
            println( "|P| = ", cardin( poly_a ) );
            println( "|QS| = ", cardin( QS ) );
            println( "|Q| = ", cardin( poly_b ) );
            println( "Computing radii simplified Frechet distance..." );
        end
        #    m_mid = frechet_ve_r_mono_compute( PS, QS  );
        f_debug && println( "Approx refinement : ", aprx_refinement );
        m_mid, f_exact, PSR, QSR = frechet_mono_via_refinement( PS, QS,
                                                              aprx_refinement );

        f_debug  &&  println( "frechet mono via refinment computed" );
        m_a = frechet_ve_r_mono_compute( poly_a, PSR );
        mmu = Morphing_combine( m_a, m_mid );
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
            continue;
        end

        # Now we compute the distance, with offsets...
        PSR_offs = Morphing_extract_offsets( m_a )[2]
        QSR_offs = Morphing_extract_offsets( m_b )[1]

        m_final = frechet_ve_r_compute_ext( PSR, QSR, PSR_offs, QSR_offs,
                                            true );
        if  ( floating_equal( m_final.leash,  mw.leash ) )
#            println( "We are done!" );
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
        if  f_accept_approx  &&  ( 1.000001 * mw.leash > m_final.leash )
            f_debug && println( "MW Leash length: ", mw.leash );
            f_debug && println( "m_final.leash  : ", m_final.leash );
#            println( "We are done!" );
            return  mw
        end
    end  # While loop end
end



#
#end

#
# End of file
##########################################################################
##########################################################################
##########################################################################

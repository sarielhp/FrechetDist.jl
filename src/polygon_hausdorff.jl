# Contains code to do hausdorff simplification of a curve

module polygon_hausdorff

using ..point;
using ..segment;
using ..polygon


"""
    hausdorff_dist_subseg

    In the following, for a polygon P, and a range i:j, let

    HD(P, i:j) denote the Hausdorff distance between the
    subpolygon specified by P[i:j] = P[i]P[i+1] … P[j], and the segment
    P[i]P[j].

    The parameter rng is a range - that is a pair of integer i:j.

    Computes HD( P, rng ).
"""
function  hausdorff_dist_subseg( P::Polygon{D,T},
                                rng::UnitRange{Int64} = 0:0
                                )::T where {D,T}
    if  ( rng == 0:0 )
        rng = 1:length( P )
    end
    if  ( length( rng ) <= 2 )  ||  ( length( P ) <= 2 )
        return  0;
    end

    s,t = P[ first( rng ) ], P[ last( rng ) ];
    leash::Float64 = 0;
    for  i  in  first(rng)+1:last(rng)-1
        leash = max( leash, dist_iseg_nn_point( s, t, P[ i ] ) )
    end

    return  leash;
end


"""
    exp_search_prefix

    Using exponential search, starting at vertex start of P, computes
    the first prefix ending at hi, (approximately) that has Hausdorff
    distance > w. Thus, the returned value is an index hi, such that

    HD( P, start:hi ) > w.

    If no such index found, hi is the index of the last vertex in P.
"""
function  exp_search_prefix( P::Polygon{D,T},
                             start::Int64,
                             w::Float64 ) where {D,T}
    len = cardin( P );
    hi::Int64 = min( start + 2, len );

    ( hi >= len )  &&  return  hi;

    while  hi < len
        r = hausdorff_dist_subseg( P, start:hi )
        ( r > w )  &&  return  min( hi + 5, len );
        hi = start + 2 * ( hi - start )
    end

    return  len;
end


"""
    h_bin_search_inner

    Given range i:j, and a *start* vertex of P, performs a binary search
    for the first index t ∈ i:j, such that HD(P[start:t]) ≤ w and
    HD(P[start:t+1]) > w. The algorithm uses binary search.
"""
function  h_bin_search_inner( P::Polygon{D,T}, start::Int64,
                              i::Int64, j::Int64, w::T ) where {D,T}

    if ( i >= j )  ||  ( ( start + 1) == j )
        return  j; # the distance is zero, nothing to do.
    end
    r = hausdorff_dist_subseg( P, start:j )
    if  ( r <= w )
        return  j;
    end
    if (i+1) == j
        return  i;
    end
    j = j - 1;
    mid = ( i + j ) >> 1;
    r_m = hausdorff_dist_subseg( P, start:mid )
    if  ( r_m > w )
        return   h_bin_search_inner( P, start,   i, mid - 1, w );
    end
    return       h_bin_search_inner( P, start, mid, j, w );
end


function  h_find_prefix( P::Polygon{D,T}, i::Int64,
    j::Int64, w::T )  where {D,T}

    r = hausdorff_dist_subseg( P, i:j )
    if  ( r <= w )
        return  j;
    end
    return  h_bin_search_inner( P, i, i, j, w );
end


function  hausdorff_simplify( P::Polygon{D,T}, w::T ) where {D,T}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    card = cardin( P );
    if  ( card == 0 )
        return pout, pindices;
    end
    push!( pout, P[1] );
    push!( pindices, 1 );

    curr_ind = 1;
    next_ind::Int64 = 1;
    while  true
        hi = exp_search_prefix( P, curr_ind, w );
        next_ind = h_find_prefix( P, curr_ind, hi, w )
        @assert( next_ind > curr_ind );
        push!( pindices, next_ind );
        push!( pout,  P[ next_ind ] );
        if  next_ind == card
            return  pout, pindices;
        end
        curr_ind = next_ind;
    end
end

export  hausdorff_simplify;
export  hausdorff_dist_subseg;

end

# Contains code to do hausdorff simplification of a curve

module polygon_hausdorff

using point;
using segment;
using ..polygon


"""
    hausdorff_dist_subseg

Computes the Hausdorff distance between the subpolygon specified by
P[rng], and they segment they define.
"""
function  hausdorff_dist_subseg( P::Polygon{D,T},
                                rng::UnitRange{Int64} = 0:0
                                )::T where {D,T}
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

    seg = Segment{D,T}( P[ first( rng ) ], P[ last( rng ) ] );

    leash::Float64 = 0;
    for  i  in  first(rng)+1:last(rng)-1
        q = nn_point( seg, P[ i ] );
        leash = max( leash, Dist( q, P[ i ] ) );
    end

    return  leash;
end


function  hausdorff_exp_search_prefix( P::Polygon{D,T},
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




function  hausdorff_prefix_inner(
    P::Polygon{D,T},
    ind_start::Int64, i::Int64, j::Int64, w::T ) where {D,T}

    if i >= j
        return  j;
    end
    if  (ind_start + 1) == j
        return  j; # the distance is zero, nothing to do.
    end
    r = hausdorff_dist_subseg( P, ind_start:j )
    if  ( r <= w )
        return  j;
    end
    if (i+1) == j
        return  i;
    end
    j = j - 1;
    mid = ( i + j ) >> 1;
    r_m = hausdorff_dist_subseg( P, ind_start:mid )
    if  ( r_m > w )
        return   hausdorff_prefix_inner( P, ind_start, i, mid - 1, w );
    end
    return  hausdorff_prefix_inner( P, ind_start, mid, j, w );
end

function  hausdorff_find_prefix( P::Polygon{D,T}, i::Int64,
    j::Int64, w::T )  where {D,T}

    r = hausdorff_dist_subseg( P, i:j )
    if  ( r <= w )
        return  j;
    end
    return  hausdorff_prefix_inner( P, i, i, j, w );
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
        hi = hausdorff_exp_search_prefix( P, curr_ind, w );
        next_ind = hausdorff_find_prefix( P, curr_ind, hi, w )
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

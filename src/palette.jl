

function   fill_plt( plt::Vector{T}, P::Polygon{N,T},
                     i::Int64, j::Int64 ) where {N,T}
    ( i > j )  &&  return 0.0;

    mid = (i + j) >> 1;
    if  ( mid == i )  ||  ( mid == j )
        return  0.0;
    end
    v_l = fill_plt( plt, P, i, mid );
    v_r = fill_plt( plt, P, mid, j );

    v = frechet_width_approx( P, i:j );

    plt[ mid ] = max( v_l, v_r );

    return  v; #max( v, plt[ mid ] );
end


function   frechet_palette( P::Polygon{N,T} )::Vector{T} where {N,T}
    crd = cardin( P );
    plt = zeros( T, crd )
    ( crd < 2 )  &&  return  plt;

    fill_plt( plt, P, 1, crd );

    return  plt;
end

function  p_extract( pout::Polygon{N,T},
                     pout_is::Vector{Int64},
                     P::Polygon{N,T}, plt::Vector{T}, i, j,
                     w::T)  where {N,T}
    ( i > j )  &&  return  0.0;
    mid = (i + j) >> 1;
    if  ( mid == i )  ||  ( mid == j )
        return 0.0;
    end
    if  plt[ mid ] < w
        push!( pout, P[ mid ] );
        push!( pout_is, mid );
        return plt[ mid ];
    end
    w_l = p_extract( pout, pout_is, P, plt, i, mid, w );
    push!( pout, P[ mid ] );
    push!( pout_is, mid );
    w_r = p_extract( pout, pout_is, P, plt, mid, j, w );
    #println( "w_l :", w_l );
    return  max( w_l, w_r );
end


"""
    frechet_palette_approx

Extract an approximation to P using the palette, and a width
parameter. The output polygon Q has the property the Frechet distance
betwwen P and Q is at most w (might be smaller).

The approximation provided is fast (i.e., running time proportional to
the nubmer of vertices of Q, but the approximation might be of low
quality.

"""
function   frechet_palette_approx( P::Polygon{N,T},
                                   plt::Vector{T},
                                   w::T
                                        ) where {N,T}
    len = length( plt );

    @assert( len == cardin( P ) );

    pout = Polygon{N,T}();
    pout_is = Vector{Int64}();
    if  ( cardin( P ) <= 4 )
        for  i in 1:len
            push!( pout, P[ i ] );
            push!( pout_is, i );
        end
        return  pout, pout_is, zero( T );
    end

    push!( pout, P[ 1 ] );
    push!( pout_is, 1 );

    w = p_extract( pout, pout_is, P, plt, 1, len, w );

    push!( pout, P[ len ] );
    push!( pout_is, len );

    return  pout, pout_is, w;
end


function  p_extract_level( pout::Polygon{N,T},
                     pout_is::Vector{Int64},
                     P::Polygon{N,T}, plt::Vector{T},
                     i::Int64, j::Int64,
                     level::Int64 )  where {N,T}
    if ( ( i > j )  ||  ( level <= 0 ))
        return;
    end
    
    mid::Int64 = (i + j) >> 1;
    if  ( mid == i )  ||  ( mid == j )
        return;
    end
    p_extract_level( pout, pout_is, P, plt, i, mid, level - 1 );
    push!( pout, P[ mid ] );
    push!( pout_is, mid );
    p_extract_level( pout, pout_is, P, plt, mid, j, level - 1 );
end


"""
    frechet_palette_level

Extract the approximate polygon by specifying a level of the
approximation. 

"""
function frechet_palette_level( P::Polygon{N,T},
    plt::Vector{T}, level::Int64 ) where {N,T}

    len = cardin( P );
    @assert( len > 0 );

    pout = Polygon{N,T}();
    pout_is = Vector{Int64}();
    
    push!( pout, P[ 1 ] );
    push!( pout_is, 1 );

    p_extract_level( pout, pout_is, P, plt, 1, len, level );

    push!( pout, P[ len ] );
    push!( pout_is, len );
    println( "|pout|: ", cardin( pout ) );
    return  pout, pout_is;
end

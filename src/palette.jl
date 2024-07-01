

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
0

function   frechet_approx_from_palette( P::Polygon{N,T},
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

using Parameters
using DataStructures
using Printf



mutable struct  PolygonHierarchy{D,T}
    len::T;
    P::Polygon{D,T};
    plt::Vector{T};

    polys::Vector{Polygon{D,T}};
    widths::Vector{Float64};
    lock::ReentrantLock;
end

function ph_push!( ph::PolygonHierarchy{D,T}, Q::Polygon{D,T}, w::Float64
                   )   where {D,T}
    @assert( cardin( Q ) > 1 )

    lock( ph.lock );
    if  ( length( ph.widths ) > 0 )
        if  ( ( last( ph.widths ) >= w )
            && ( cardin( last( ph.polys ) ) >= cardin(Q) ) )
            pop!( ph.polys );
            pop!( ph.widths );
        end
    end
    Ps = ph.polys;
    Ws = ph.widths
    push!( Ps, Q );
    push!( Ws, w );

    # sort in decreasing order...
    i = length( ph.polys )
    while  ( i > 1 )  &&   ( Ws[ i - 1 ] < Ws[ i ] )
        Ws[ i - 1 ], Ws[ i ] = Ws[ i ], Ws[ i - 1 ]
        Ps[ i - 1 ], Ps[ i ] = Ps[ i ], Ps[ i - 1 ]
        i = i - 1;
    end
    unlock( ph.lock );
end


function  ph_push_target( ph::PolygonHierarchy{D,T}, w::Float64,
    P::Polygon{D,T}, lmt::Int64, w_extra::Float64 )  where {D,T}

    R, R_indices = frechet_simplify_w_exp( P, w );
    if  ( cardin( R ) > lmt )
        return  -1.0;
    end

    mR = frechet_c_mono_approx_subcurve( P, R, R_indices )[ 1 ];
    new_w = mR.leash + w_extra;

    ph_push!( ph, R, new_w );

    return  new_w;
end


function  ph_print( ph::PolygonHierarchy{D,T} ) where  {D,T}
    P = ph.P
    println( "--------------------------------------------" );
    println( "|P|: ", cardin( P ) );
    for  i  in  1:length( ph.widths )
        println( "w[", i, "] :", ph.widths[ i ],
                 "  #: ", cardin( ph.polys[ i ] ) );
    end
    println( "----------------------------------------------\n" );
end


function   ph_find_linear( Ws::Vector{Float64}, w, resolution )
    for  i  in  1:length( Ws )
        #println( ph.widths[ i ] );
        if  ( Ws[ i ] <= w )
            if  ( w <= ( resolution * Ws[ i ] ) )
                #println( "Bingo!" );
                return  i
            else
                return  -1;
            end
        end
    end
    return  -1;
end
function   ph_approx_binary_search( Ws, i, j, w, resolution )
    while  true
        #println( "[",i,"...",j, "]" );
        # w[i] >= w >= w[j]
        ( i > j )       &&  return  -1;
        if  ( Ws[ i ] <= w)
            if  ( w <= Ws[ i ]*resolution )
                return  i;
            else
                return  -1;
            end
        end

        # Ws[i] > w
        if  ( i == j )
            return  -1;
        end

        mid = ( i + j ) >> 1;
        if  Ws[ mid ] >= w
            i = mid + 1;  # continuing into mid+1:j
            continue;
        end
        # Ws[mid ] < w
        j = mid;
    end
end

"""
    ph_approx

Tries to find an approximation of the curve with width at most w that
was already computed, as long as w is not smaller than
w/resolution. If such a simplification was found, then we try to find
the best simplification with at most n_max vertices.

Otherwise, it computes such a simplification, stores it in the
hierarchy, and returns it.

"""
function   ph_approx( ph::PolygonHierarchy{D,T}, w::Float64,
                      n_max::Int64 = 25,
                      resolution::Float64 = 4.0 ) where  {D,T}
    f_verify::Bool = false

    i = ph_find_linear( ph.widths, w, resolution );
    if  ( i > 0 )
        len = length( ph.polys );
        while  ( i < len )  &&  ( cardin( ph.polys[ i + 1 ] ) < n_max )
            i = i + 1;
        end
        return  ph.polys[ i ], ph.widths[ i ]
    end

    new_w = w / 4.0;
    i = ph_find_linear( ph.widths, new_w, resolution );
    if  ( i > 0 )
        Z, wZ = ph.polys[ i ],  ph.widths[ i ]
    else
        Z, grb, wZ = frechet_palette_approx( ph.P, ph.plt, new_w );
    end

    X, X_indices = frechet_simplify_w_exp( Z, max( w-wZ, 0.0 ) );
    mX = frechet_c_mono_approx_subcurve( Z, X, X_indices )[ 1 ];

    w_out = (mX.leash + wZ);

    #m = frechet_c_compute( X, ph.P );
    #println( cardin( ph.P ), " ",cardin( X ), " err: ", w_out, " rerr: ", m.leash );
    #=
    if  ( ( cardin( X ) == 28 )  &&  (cardin( ph.P ) > 200 ) )
        Polygon_write_to_file( X, "X.txt" );
        Polygon_write_to_file( ph.P, "P.txt" );
        exit( -1 );
    end
    =#
    ph_push!( ph, X, w_out );

    return  X, w_out;
end


function  ph_compute_hierarchy( ph::PolygonHierarchy{D,T},
                                ratio::Float64 = 4.0,
                                lmt::Int64 = 500 )  where  {D,T}
    P = ph.P;
    lmt::Int64 = min( round( Int64, cardin( P ) * 0.95 ), lmt );

    w = ph.widths[ 1 ];
    for  i in 1:200
        if cardin( last( ph.polys ) ) >= lmt
            break
        end
        wA = w / ratio;

        Z, Z_indices, wZ = frechet_palette_approx( P, ph.plt, wA / 4.0 );

        w = ph_push_target( ph,   max( wA - wZ, 0.0 ), Z, lmt, wZ )
        ( w < -0.5 )  &&  break;
    end
end


function ph_init( P::Polygon{D,T}, ph_lock ) where  {D,T}
    ph = PolygonHierarchy{D,T}( 0.0, P, Vector{Int64}(), Vector{Polygon2F}(),
                           Vector{Float64}(), ph_lock );
    ph.len = Polygon_length( P );
    ph.plt = frechet_palette( P );

    w = frechet_width_approx( P );
    ph_push!( ph, Polygon_spine( P ), w );

    return  ph;
end

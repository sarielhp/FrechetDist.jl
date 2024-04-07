# Originally contributed by S. Har-Peled
# under MIT License

#---------------------------------------------------------------
# Code to compute discrete Frechet. Both the standard type,
# and the retractable type. Note, that the retractable version
# still allocates the quadratic size table (unlike the ve version
# which uses hashing to avoid this).
#---------------------------------------------------------------


function   d_frechet_extract_solution( P::Polygon{D,T}, Q,
    dp_dec_i, n_p, n_q )::Morphing{D,T} where {D,T}
    i = n_p;
    j = n_q;

    #println( "Extract solution?" );
    up::Int64 = n_p + n_q + 1;
    peout = Vector{EventPoint{D,T}}( undef, up );
    qeout = Vector{EventPoint{D,T}}( undef, up );

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



function   frechet_d_compute_inner( P::Polygon{D,T}, Q::Polygon{D,T},
                                    dp::Array{Float64, 2}, dp_dec_i
                                        ) where {D,T}
    d::Float64 = 0;
    jp::Int64 = 0;
    ip::Int64 = 0;

    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    iters::Int64 = 0;
    dp[ 1, 1 ] = Dist( P[ 1 ], Q[ 1 ] );
    for  i::Int64 in 1:n_p
        ip = max( i - 1, 1 );
        for  j in 1:n_q
            iters = iters + 1;
            d = Dist( P[ i ], Q[ j ] );
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

    return  iters;
end


"""
    frechet_d_compute

Computes the discrete Fréchet distance between the two sequences of
points. It interpresents the vertices of the polygons and the desired
sequence of points.

"""
function   frechet_d_compute( P::Polygon{D,T},
                              Q::Polygon{D,T} ) where {D,T}
    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    dp::Array{Float64, 2} = Array{Float64, 2}(undef, n_p,n_q);
    dp_dec_i = falses( n_p, n_q );

    iters = frechet_d_compute_inner( P, Q, dp, dp_dec_i );

    m = d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
    m.iters = iters;

    return  m;
end


function   frechet_d_compute_dist( P::Polygon{D,T},
                                   Q::Polygon{D,T} ) where {D,T}
    n_p::Int64 = cardin( P );
    n_q::Int64 = cardin( Q );

    dp::Array{Float64, 2} = Array{Float64, 2}(undef, n_p,n_q);
    dp_dec_i = falses( n_p, n_q );

    frechet_d_compute_inner( P, Q, dp, dp_dec_i );

    return  dp[ n_p, n_q ];
end



##########################################################################
##########################################################################
# _lopt_frechet
"""
    frechet_d_r_compute

Compute discrete frechet distance that is locally optimal
Frechet. Essentially discrete frechet + Prim/Dijkstra algorithm For
the discrete case, that is modified to be retractable -- that is
minimize the maximum bottleneck edge being computed.
"""
function   frechet_d_r_compute( P::Polygon{D,T}, Q::Polygon{D,T}
                                            ) where {D,T}
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

    iters::Int64 = 0;
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

    m = d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
    m.iters = iters;

    return  m;
end


function   DTW_d_compute( P::Polygon{D,T}, Q::Polygon{D,T}
                          ) where {D,T}

    pq = PriorityQueue{Tuple{Int64, Int64},Float64}();

    n_p = cardin( P );
    n_q = cardin( Q );

    dp = zeros( n_p, n_q );
    dp_dec_i = falses( n_p, n_q );
    handled = falses( n_p, n_q );

    function  push_value( ia, ja, val, i )
        new_val = val + Dist( P[ ia ], Q[ ja ] )

        if  ( ( dp[ ia, ja ] > 0.0 )
              &&  ( dp[ ia, ja ] <= new_val ) )
            return;
        end
        dp[ ia, ja ] = new_val;
        pq[ ia, ja ] = new_val;  # Push to queue
        dp_dec_i[ ia, ja ] = ( i < ia );
    end

    dp[ 1, 1 ] = Dist( P[ 1 ] , Q[ 1 ] );
    #handled[ 1, 1 ] = true;

    enqueue!( pq, (1,1), dp[ 1, 1 ] );

    iters::Int64 = 0;
    while  ! isempty( pq )
        iters = iters + 1;

        ele = peek( pq );
        i, j = ele[ 1 ];
        value = ele[ 2 ];

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
            push_value( i + 1, j, value, i );
        end
        if  ( j < n_q )
            push_value( i, j+1, value, i );
        end
    end

    m = d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
    m.iters = iters;

    #exit( -1 );
    return  m;
end



###################################################################

"""
    ADTW_lb_compute

Compute a lower bound on the ADTW distance between two curves P and Q.

"""
function   ADTW_lb_compute( P::Polygon{D,T}, Q::Polygon{D,T}
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

        if  j > 1
            lb = min( lb,
                iseg_iseg_dist( PA[i ], PA[i + 1], PB[ j - 1 ], PB[ j ] )
            );
        end
        if  j < cardin( PB )
            lb = min( lb,
                iseg_iseg_dist( PA[ i ], PA[i + 1], PB[ j ], PB[ j + 1 ] )
            );
        end
        return  lb;
    end

    function edge_v( PA::Polygon{D,T}, i::Int64, PB::Polygon{D,T},
        j::Int64 )::Float64 where {D,T}
        if  ( i == cardin( PA ) )
            return 0.0;
        end
        #println( "edge_lb( ?, ", i, ", ?, ", j, ")" ); 
        lb = edge_lb( PA, i, PB, j );
        return  Dist( PA[ i ], PA[ i + 1 ] ) * lb;
    end


    function  push_value( ia::Int64, ja, new_val::Float64, i )
        #new_val = val + Dist( P[ ia ], Q[ ja ] )

        if  ( ( dp[ ia, ja ] > 0.0 )
              &&  ( dp[ ia, ja ] <= new_val ) )
            return;
        end
        dp[ ia, ja ] = new_val;
        pq[ ia, ja ] = new_val;  # Push to queue
        dp_dec_i[ ia, ja ] = ( i < ia );
    end

    dp[ 1, 1 ] = 0.0;
    #handled[ 1, 1 ] = true;

    enqueue!( pq, (1,1), dp[ 1, 1 ] );

    iters::Int64 = 0;
    while  ! isempty( pq )
        iters = iters + 1;

        ele = peek( pq );
        i::Int64, j::Int64 = ele[ 1 ];
        value = ele[ 2 ];

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
            push_value( i + 1, j, value + edge_v( P, i, Q, j ), i );
        end
        if  ( j < n_q )
            push_value( i, j+1, value + edge_v( Q, j, P, i ), i );
        end
    end

    m = d_frechet_extract_solution( P, Q, dp_dec_i, n_p, n_q );
    m.iters = iters;

    m.sol_value = dp[ n_p, n_q ];

    #exit( -1 );
    return  m;
end




"""
    frechet_d_compute_sample

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

- f_lopt
     true: Use the rectractable Frechet distance version
     false: Use the standard discrete Frechet version.
"""
function   frechet_d_r_compute_sample( polya::Polygon{D,T},
                                       polyb::Polygon{D,T},
                                       n::Int64,
                                       f_lopt::Bool = true
                                     )    where {D,T}
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


#
# End of file
##########################################################################

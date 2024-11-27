module segment

using StaticArrays
using Parameters
using LinearAlgebra
using DelimitedFiles

using point;


###############################################33
###############################################33
### Line type
###############################################33


"""
    Line in D dimensions.

    `p` is a point on the line and `u` is the direction vector (not
    necessarily normalized). Parametrised as \$p + ut\$

"""
struct Line{D,T}
    p::MVector{D,T}
    u::MVector{D,T}
end

###############################################33
###############################################33
### Segment type
###############################################33
#Point{D,T} = point.Point{D,T}

"""
    Segment

Specifies a *directed* segment by two endpoints.
"""
struct Segment{D,T}
    p::Point{D,T}
    q::Point{D,T}
end

function  Segment{D,T}() where{D,T}
    return   Segment{D,T}(Point{D,T}(), Point{D,T}() );
end

@inline function  at( seg::Segment{D,T}, t::Float64 ) where{D,T}
    return convex_comb( seg.p, seg.q, t );
end



function Base.show(io::IO, s::Segment{D,T}) where {D,T}
    print( io," [[" );
    print( io, s.p );
    print( io, ".." )
    print( io, s.q );
    print( io, "]] " );
end


@inline function  Segment_length( seg::Segment{D,T} ) where {D,T}
    return point.Dist( seg.p, seg.q );
end

@inline function  Dist( s::Segment{D,T}, qr::Point{D,T} ) where {D,T}
    nn = nn_point( s, qr );
    return  point.Dist( nn, qr );
end

@inline function  Dist( a::Segment{D,T}, b::Segment{D,T} ) where {D,T}
    return  iseg_iseg_dist( a.p, a.q, b.p, b.q );
end

@inline function  Segment_get_convex_coef( s::Segment{D,T}, qr::Point{D,T} ) where{D,T}
    d = point.Dist( s.p, qr );
    len = Segment_length( s );
    if  len == 0.0
        return  0.0;
    end
    return  d/len;
end

"""
    iseg_nn_point_ext_ext

    Two parameters should be precomputed before calling this function:

   seg_len_sq = DistSq( s_p, s_q ) 
   s_vec =  sub(s_q, s_p)

   return  nearest point, and its convex combination coefficient;
"""
@inline function  iseg_nn_point_ext_ext( s_p::Point{D,T}, s_q::Point{D,T},
    qr::Point{D,T}, seg_len_sq::T, s_vec::Point{D,T} ) where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    # D(t) = ||u(t)||^2 = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    #
    # s_vec = sub(s_q, s_p)
    a = seg_len_sq #DistSq( s_p, s_q );
    if   ( a == 0.0 )
        return  s_p, zero(T);
    end
    b = 2.0 * point.dot( sub(s_p, qr), s_vec );
    t::T = -b / (2.0 * a);

    if  ( t <= 0.0 )
        return  s_p, zero(T);
    elseif  ( t >= 1.0 )
        return  s_q, one(T);
    end

    return  convex_comb( s_p, s_q, t ), t;
end

# Get nearest point on segment...
"""
    iseg_nn_point_ext

Returns the closest point to the segment induced by the first two
points, to the query point. By avoiding creating the segment iself, it
is hopefully more efficient. Returns the convex_combination param of the nn
point.
"""
@inline function  iseg_nn_point_ext( s_p::Point{D,T}, s_q::Point{D,T},
                                qr::Point{D,T} ) where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    # D(t) = ||u(t)||^2 = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    a = DistSq( s_p, s_q );
    if   ( a == 0.0 )
        return  s_p, zero(T);
    end
    b = 2.0 * point.dot( sub(s_p, qr), sub(s_q, s_p) );
    t::T = -b /(2.0 * a);

    if  ( t <= 0.0 )
        return  s_p, zero(T);
    elseif  ( t >= 1.0 )
        t = 1;
        return  s_q, one(T);
    end
    pon = convex_comb( s_p, s_q, t );

    return  pon, t;
end


function  iseg_nn_point( s_p::Point{D,T}, s_q::Point{D,T},
                                qr::Point{D,T} ) where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    # D(t) = ||u(t)||^2 = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    a = DistSq( s_p, s_q );
    if   ( a == 0.0 )
        return  s_p;
    end
    b = 2.0* point.dot( sub(s_p, qr), sub(s_q, s_p) );
    t::T = -b /(2.0 * a);

    if  ( t <= 0.0 )
        return  s_p;
    end
    if  ( t >= 1.0 )
        return  s_q;
    end
    pon = convex_comb( s_p, s_q, t );

    return  pon;
end


function  iseg_iseg_dist( a_p::Point{D,T}, a_q::Point{D,T},
    b_p::Point{D,T}, b_q::Point{D,T} ) where {D,T}
    # v(s,t) = a_p*(1-s) + a_q * s - b_p*(1-t) - b_q * t
    #      = (a_p - b_p)  +  s * (a_q-a_p)  +  t * (b_p-b_q)
    #      = v_1  +  s * v_2  +  t * v_3
    #v_1::Point{D,T};
    #v_2::Point{D,T};
    #v_3::Point{D,T};

    v_1 = sub(a_p, b_p);
    v_2 = sub(a_q, a_p);
    v_3 = sub(b_p, b_q);

    #=

    D(s,t)
    = <v_1 + s*v_2 + t *v_3, v_1 + s*v_2 + t *v_3>
    = ||v_1||^2  +  2 * s * <v_1, v_2>  +  2 * t *<v_1,v_3>
       + s^2 * ||v_2||^2 + 2*s*t*<v_2,v_3>  + t^2 ||v_3||^2
    =
    Need to solve the linear system:

    0 = D'_s(s,t) = 2*<v_1,v_2> + 2*s*||V_2||^2  +  2*t<v_2,v_3>
    0 = D'_t(s,t) = 2*<v_1,v_3> + 2*s*<v_2,v_3> + 2*t * ||v_3||^2

    or equivalently:

    -<v_1,v_2> = s * ||V_2||^2  +  t * <v_2,v_3>
    -<v_1,v_3> = s * <v_2,v_3>  +  t * ||v_3||^2

    =#

    # Coefficients
    c = [ (- point.dot( v_1, v_2 )),   (-point.dot( v_1, v_3 ) ) ];

    m = [ ( point.dot(v_2, v_2))  ( point.dot(v_2, v_3) );
          ( point.dot(v_2, v_3))  ( point.dot(v_3, v_3) ) ];

    rk = rank( m );
    @assert( rk > 0 );
    s::Float64 = 0.0;
    y::Float64 = 0.0;
    #println( "rk = ", rk );
    if  ( rk == 1 )
        # The minimum distance is realized by one of the endpoints.
        d = min(
            dist_iseg_nn_point( a_p, a_q, b_p ),
            dist_iseg_nn_point( a_p, a_q, b_q ),
            dist_iseg_nn_point( b_p, b_q, a_p ),
            dist_iseg_nn_point( b_p, b_q, a_q ) );
        return  d;
    end
    #=
        println( "--- a, b --------------------" );
        println( a_p, a_q );
        println( b_p, b_q );

        println( "--- v_1,2,3 --------------------" );

        println( v_1 );
        println( v_2 );
        println( v_3 );

        println( "\n\n--- m --------------------" );
        println( "rank m: ", rank( m ) );
        println( m );
        println( "\n\n--- c --------------------" );
        println( c );

        println( "m.size: ", size( m ) );
        println( "c.size: ", size( c ) );
        =#

    # Solve the system...
    b = m \ c;
    #println( "b.size: ", size( b ) );
    s = b[ 1 ];
    t = b[ 2 ];

    # Snap solution if needed to the [0,1.0] interval....
    s = max( min( s, 1.0 ), 0.0 )
    t = max( min( t, 1.0 ), 0.0 )

    d::Float64 = point.Dist( convex_comb( a_p, a_q, s ),
        convex_comb( b_p, b_q, t ) );

    d = min( d, dist_iseg_nn_point( a_p, a_q, b_p ),
                dist_iseg_nn_point( a_p, a_q, b_q ),
                dist_iseg_nn_point( b_p, b_q, a_p ),
                dist_iseg_nn_point( b_p, b_q, a_q ) );


    da = point.Dist( a_p, b_p );
    if  ( da < d )  &&  ( abs( da - d ) > (1e-8 * (d + da ) ) )
        println( "da < d? " );
        println( "da: ", da );
        println( "d ", d );
        println( "s: ", s );
        println( "t: ", t );
        println( "b: ", b );
        exit( -1 );

    end
    db = point.Dist( a_q, b_q );
    if  ( db < d )  &&  ( abs( db - d ) > (1e-8 * (d + db ) ) )
        println( "db < d? " );
        println( "db: ", da );
        println( "d ", d );
    end

    #println( "d = ", d );

    return  d;
end



#############################################
# Get nearest point on segment...
#
# s_p - s_q : the segment
# qr : the query point.
#############################################
"""
    dist_iseg_nn_point

Returns the *distance* to the closest point lying on the segment induced by
the first two points, to the query point. By avoiding creating the
segment iself, it is hopeflly more efficient.
"""
@inline function dist_iseg_nn_point( s_p::Point{D,T},
                                     s_q::Point{D,T}, qr::Point{D,T}
)  where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    #
    # D(t) = ||u(t)||^2
    #      = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    #
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    #dff = Point{D,T}()

    #sub_dst( s_p, qr, dff );
    dff = sub( s_p, qr );

    a = DistSq( s_p, s_q );
    b = 2.0* point.dot( dff, sub( s_q, s_p ) );

    c = point.dot( dff, dff );

    t = -b /(2.0 * a);

    local val::Float64;
    if  ( t < 0 )
        #t = 0;
        val = c ;
    elseif  ( t > 1 )
        #t = 1;
        val = a + b + c;
    else
        val = a*t^2 + b*t + c;
    end
    if  ( val <= 1e-20 )
        return 0;
    end
    return sqrt( val );
end


"""
    nn_point

Returns the closest point on the segment `s` to the query point `qr`.

"""
@inline function  nn_point( s::Segment{D,T}, qr::Point{D,T} ) where {D,T}
    return   iseg_nn_point( s.p, s.q, qr );
end


#######################################################################
# Check if the plane bisector of p and q intersect seg, and if so where...
#f_on,t,p
#######################################################################
"""
    Segment_get_bisection_point -> Bool, Real, Point

    Computes the intersection point of the segment `seg` with the
bisector plane between `p` and `q`.

# Returns

The first argument returns whether the segment intersects the
bisector, the pramaterized location (tm), and the intersection piont
itself.

"""
function  Segment_get_bisection_point( seg::Segment{D,T}, p::Point{D,T},
    q::Point{D,T}  ) where {D,T}
    # Consider the segment going through p and q, and its middle point (mid).
    dir = sub( q, p );
    mid = ( q + p ) /2;
    pos = point.dot( dir, mid );

    # The line induced by seg, we are interested in the t, when
    # its dot proct with dir is equal to pos.

    # p::Point{D,T}
    # q::Point{D,T}
    vec = seg.q - seg.p;

    # seg.p + tm * vec:  The line under consideration
    # So:
    #     dot( dir, seg.p ) + t * dot( dir, vec ) = pos
    tm::Float64 = ( pos - point.dot( dir, seg.p ) ) / point.dot( dir, vec ) ;

    f_on::Bool = (0.0 <= tm <= 1.0);
    out = at( seg, tm );

#    println( "" );
#    println( "DDDD 1:", Dist( p, out ), "  D2: ", Dist( q, out ) );
#    println( "" );

    return  f_on, tm, out;
end

Segment2F = Segment{2,Float64};


export  Line;
export  Segment, Segment2F;

export  Segment_get_bisection_point
export  nn_point
export  at, Segment_get_convex_coef
export  Segment_length
export  Segment_get_bisection_point
export  dist_iseg_nn_point
export  iseg_nn_point, iseg_nn_point_ext, iseg_nn_point_ext_ext
export  iseg_iseg_dist

end # module segment

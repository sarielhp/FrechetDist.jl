#
# DistFunction is (informally) a function encoding the distance of a
# point moving on a line, to its nearest point on another line.
#•
# More formally, it is a function of the form:
#
# f(x) = sqrt( a*x^2 + b*x + c ) where a,b,c are constants, and the
# function is well defined for all x.
#
#
using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles


mutable struct DistFunction{T}
    a::T;
    b::T;
    c::T;
end

function  eval( f::DistFunction{T}, x::T ) where {T}
    return  sqrt( f.a*x*x + f.b*x + f.c );
end

function  integral( f::DistFunction{T}, x::T ) where {T}
    v = eval( f, x );

    t = (2.0*f.a*x + f.b ) / sqrt(4.0*f.a*f.c-f.b*f.b);

    as = asinh( t );

    A = f.c * as / ( 2.0*sqrt( f.a ) );
    B = - (  f.b^2 / ( 8.0*f.a^(3/2) )  ) * as;
    C = (x/2.0) * v;
    D = (f.b/(4.0 * f.a) ) * v;

    return A + B + C + D;
end


###################################################################3
# Since I got the integral of the function f from software (i.e.,
# sagemath), and I was too lazy to verify it by hand, I did a
# numerical verification.
###################################################################3

function  DistFunc_verify_integral( a, b, c, x )
    f = DistFunction( a, b, c );

    delta = x/2.0;
    local v, diff;
    for i in 1:14
        l = x - delta;
        hf = x + delta;
        diff = ( integral( f, h ) - integral( f, l )  ) / (2.0*delta);
        v = eval( f, x )
        delta = delta / 2.0;
    end
    print( "f (", x,"): ", v, "        " );
#    println( "F'(", x,"): ", diff  );
    println( "f - F' ",v - diff );
end

function DistFunc_verify_integral_2( a, b, c )
    DistFunc_verify_integral( a, b, c, 0.354 );
    DistFunc_verify_integral( a, b, c, 0.1354 );
    DistFunc_verify_integral( a, b, c, 1.1354 );
    DistFunc_verify_integral( a, b, c, 5.1354 );
    DistFunc_verify_integral( a, b, c, 7.1354 );
    DistFunc_verify_integral( a, b, c, 17.254 );
    DistFunc_verify_integral( a, b, c, 15.254 );
end

function  DistFunc_verify_integral_ext()
    DistFunc_verify_integral_2( 1.0, 1.0, 1.0 );

    DistFunc_verify_integral_2( 1.0, 1.0, 10.0 );
    DistFunc_verify_integral_2( 1.0, 1.0, 100.0 );
    DistFunc_verify_integral_2( 1.0, 2.0, 100.0 );
end


function  SweepDist_segs_p( p_a::Point{D,T},
                      p_b::Point{D,T},
                      q_a::Point{D,T},
                      q_b::Point{D,T} ) where {D,T}

    len_edge_p = Dist( p_a, p_b );

    if  ( len_edge_p == 0.0 )
        return  0.0;
    end

    @assert( len_edge_p > 0.0 );
    v_start::Point{D,T} = q_a - p_a;
    v_end = q_a - p_a;
    diff::Point{D,T} =  v_end - v_start;

    p_diff::Point{D,T} = diff / len_edge_p;
    # Location of the point at time t along the edge p.
    # p(t) := (t / len_edge_p) * diff + v_start;
    # p(t) := t  * p_diff + v_start;
    #
    # f(t) := ||p(t)|| = <p(t), p(t)>
    #       = t^2 * ||p_diff||^2 + 2*<p_diff,v_start> + <v_start,v_start>
    a = cg.dot( p_diff, p_diff );
    b = 2.0 * cg.dot( p_diff, v_start );
    c = cg.dot( v_start, v_start );

    f = DistFunction( a, b, c );

    return  integral( f, len_edge_p ) -  integral( f, 0.0 )
end

function  SweepDist_segs( p_a::Point{D,T},
                      p_b::Point{D,T},
                      q_a::Point{D,T},
                      q_b::Point{D,T} ) where {D,T}
    lp = SweepDist_segs_p( p_a, p_b, q_a, q_b );
    lq = SweepDist_segs_p( q_a, q_b, p_a, p_b );

    return  lp + lq;
end

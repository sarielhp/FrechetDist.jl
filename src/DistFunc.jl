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


#@with_kw
mutable struct DistFunction{T}
    a::T;
    b::T;
    c::T;

    disc::T;  # discriminant
    f_linear::Bool;
    root::T;
end

function DFunction( _a::T = 0.0, _b::T = 0.0, _c::T = 0.0,
                       _disc::T = 0.0, _f_linear::Bool = false,
                       _root::T = 0.0
                       ) where {T}
    return  DistFunction( _a, _b, _c, _disc, _f_linear, _root );
end

function normalize( f::DistFunction{T} ) where {T}
    f.disc = f.b*f.b - 4.0*f.a*f.c;
    if  ( abs( f.disc ) < 1e-13 )
        f.f_linear = true;
        root = -f.b/(2.0*f.a);
    end
    return  f;
end

function  dist_func( _a, _b, _c )
    f = DFunction( _a, _b, _c );

    return  normalize( f );
end




function  eval( f::DistFunction{T}, x::T ) where {T}
    v = f.a*x*x + f.b*x + f.c;
    if  ( abs( v ) < 1e-16 )
        return  0;
    end
    return  sqrt( v );
end

function  integral( f::DistFunction{T}, x::T ) where {T}
    f_debug::Bool = false;
    v = eval( f, x );

    f_debug && println( "FFF f_a: ", f.a );
    f_debug && println( "FFF f_b: ", f.b );
    f_debug && println( "FFF f_c: ", f.c );
    tmp_a = 4.0 * f.a * f.c;
    tmp_b = f.b*f.b;

    val::Float64 = tmp_a - tmp_b;
    f_debug  &&  println( "tmp_a: ", tmp_a );
    f_debug  &&  println( "tmp_b: ", tmp_b );
    f_debug  &&  println( "__VVVVVVV: ", val );
    if  ( abs( val ) < 1e-10 )
        val = 0.0;
    end

    f_debug  &&  println( "VVVVVVV: ", val );
    local as;
    if  ( val == 0.0 )
        as = 64;
    else
        t = (2.0*f.a*x + f.b ) / sqrt( val );
        as = asinh( t );
    end
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

function  integral_interval( f::DistFunction{T}, low, high ) where {T}
    @assert( low <= high );

    if  ( low == high )
        return  0;
    end
    if  ( f.f_linear )
        if  ( ( f.root <= low )  ||  ( f.root >= high ) )
            val = (high - low) * ( eval( f, low )
                                   + eval( f, high ) ) / 2.0;
            return  val;
        end
        v_l = ( f.root - low ) * eval( f, low ) / 2.0;
        v_r = ( high - f.root ) * eval( f, high ) / 2.0;
        return  v_l + v_r;
    end
    return  integral( f, high ) - integral( f, low );
    #ymin::Float64 = integral( f, 0.0 );

end


function  SweepDist_segs_p( p_a::Point{D,T},
                      p_b::Point{D,T},
                      q_a::Point{D,T},
                      q_b::Point{D,T} ) where {D,T}

    f_debug::Bool = false;
    len_edge_p = Dist( p_a, p_b );

    if  ( len_edge_p == 0.0 )
        return  0.0;
    end

    f_debug  &&  println( "len: ", len_edge_p );
    @assert( len_edge_p > 0.0 );

    f_debug  &&  println( "q_a: ", q_a );
    f_debug  &&  println( "q_b: ", q_b );
    f_debug  &&  println( "p_a: ", p_a );
    f_debug  &&  println( "p_b: ", p_b );
    v_start = q_a - p_a ;
    v_end = q_b - p_b;

    if  Dist( v_start, v_end ) < 0.0000001 * cg.norm( v_start )
        return  len_edge_p * cg.norm( v_start );
    end

    f_debug  &&  println( "v_start: ", v_start );
    f_debug  &&  println( "v_end: ", v_end );
    diff::Point{D,T} =  v_end - v_start;

    f_debug  &&  println( "Diff: ", diff );
    p_diff::Point{D,T} = (1.0 / len_edge_p) * diff;
    f_debug  &&  println( "p_diff: ", p_diff );
    # Location of the point at time t along the edge p.
    # p(t) := (t / len_edge_p) * diff + v_start;
    # p(t) := t  * p_diff + v_start;
    #
    # f(t) := ||p(t)|| = <p(t), p(t)>
    #       = t^2 * ||p_diff||^2 + 2*<p_diff,v_start> + <v_start,v_start>
    a = cg.dot( p_diff, p_diff );
    b = 2.0 * cg.dot( p_diff, v_start );
    c = cg.dot( v_start, v_start );

    f_debug  &&  println( "a: ",  a, " b: ", b, " c: ", c );
    f = dist_func( a, b, c );

    vl = integral_interval( f, 0.0, len_edge_p );

    #ymax::Float64 = integral( f, len_edge_p );
    #ymin::Float64 = integral( f, 0.0 );
    #vl::Float64 = ymax - ymin;
    #println( "ymin: ", ymin );
    #println( "ymax: ", ymax );

    if  ( isnan( vl ) )
        @assert( ! isnan( vl ) );
        println( "VL: ", vl );
    end
    f_debug  &&  println( "VL: ", vl );

    return  vl;
end

function  SweepDist_segs( p_a::Point{D,T},
                      p_b::Point{D,T},
                      q_a::Point{D,T},
                      q_b::Point{D,T} ) where {D,T}
    lp = SweepDist_segs_p( p_a, p_b, q_a, q_b );
    lq = SweepDist_segs_p( q_a, q_b, p_a, p_b );

    return  lp + lq;
end

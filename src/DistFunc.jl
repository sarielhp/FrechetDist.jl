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

function  verify( a, b, c, x )
    f = DistFunction( a, b, c );

    delta = x/2.0;
    local v, diff;
    for i in 1:14
        l = x - delta;
        h = x + delta;
        diff = ( integral( f, h ) - integral( f, l )  ) / (2.0*delta);
        v = eval( f, x )
        delta = delta / 2.0;
    end
    print( "f (", x,"): ", v, "        " );
#    println( "F'(", x,"): ", diff  );
    println( "f - F' ",v - diff );
end

function verify_several( a, b, c )
    verify( a, b, c, 0.354 );
    verify( a, b, c, 0.1354 );
    verify( a, b, c, 1.1354 );
    verify( a, b, c, 5.1354 );
    verify( a, b, c, 7.1354 );
    verify( a, b, c, 17.254 );
    verify( a, b, c, 15.254 );
end

verify_several( 1.0, 1.0, 1.0 );

verify_several( 1.0, 1.0, 10.0 );
verify_several( 1.0, 1.0, 100.0 );
verify_several( 1.0, 2.0, 100.0 );


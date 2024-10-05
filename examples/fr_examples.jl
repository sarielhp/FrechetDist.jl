# Originally contributed by S. Har-Peled
# under MIT License


############################################################
# Various examples of curves for computing Frechet distance
############################################################
function  example_1()
    polya = Polygon2F( );
    polyb = Polygon2F( );

    xs = 4;
    ys = 1.5;

    push!( polya, npoint( 0.0, 0 ), npoint( xs*1.0, 0.0 ) );
    push!( polyb,
           npoint( 0.0, ys*0.1 ), npoint( xs * 0.45, ys * 0.1 ),
           npoint( xs*0.5, ys*0.9 ), npoint( xs*0.55  , ys*0.1),
           npoint( xs*1.0, ys*0.1 ) );
    return  polya, polyb;
end

############################################################
# Various examples of curves for computing Frechet distance
############################################################
function  example_32()
    polya = Polygon2F( );
    polyb = Polygon2F( );

    xs = 4;
    ys = 1.5;

    push!( polya, npoint( 0.0, 0 ), npoint( xs*1.0, 0.0 ) );
    push!( polyb,
           npoint( 0.0, ys*0.1 ), npoint( xs * 0.15, ys * 0.1 ),
           npoint( xs*0.2, ys*0.9 ), npoint( xs*0.25  , ys*0.1),
           npoint( xs*1.0, ys*0.1 ) );
    return  polya, polyb;
end

function  example_30()
    polya = Polygon2F( );
    polyb = Polygon2F( );

    xs = Float64( 4.0 );
    ys = Float64( 0.2 );

    push!( polya, npoint( 0.0, 0 ), npoint( xs*1.0, 0.0 ) );
    push!( polyb,
           npoint( 0.0, ys ), npoint( xs , ys ),
           npoint( 0.0, 2.0*ys ), npoint( xs , 2.0*ys ) );
    return  polya, polyb;
end

function  example_2()
    polya = Polygon2F( );
    polyb = Polygon2F();

    push!( polya, npoint( 0.0, 0 ), npoint( 2.0, 0.0 ) );
    push!( polyb, npoint( 0.0, 1.0),  npoint( 0.4, 0.2 ),
           npoint( 0.6, 0.8 ),
           npoint( 1.0,0.3 ), npoint( 2.0, 0.3)  );
    return  polya, polyb;
end

function  example_3()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    push!( polya, npoint( 0.0, 0 ),
           npoint( 1.1, 0.0 ),
           npoint( 1.1, 0.1 ),
           npoint( 1.0, 0.1 ),
           npoint( 1.0, 0.2 ),
           npoint( 1.2, 0.2 ),
           npoint( 1.2, 0.0 ),
           npoint( 2.0, 0.0 ) );
    push!( polyb,
           npoint( 0.0, 0.3),
           npoint( 0.4, 0.3 ),
           npoint( 0.4, 0.6 ),
           npoint( 0.3, 0.6 ),
           npoint( 0.3, 0.7 ),
           npoint( 0.5, 0.7 ),
           npoint( 0.5, 0.3 ),
           npoint( 2.0, 0.3)  );
    return  polya, polyb;
end


function  example_4()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( 0.0, 1 )
           );
    push!( polyb,
           npoint( 0.1, 0.0 )
           ,
           npoint( 1.1, 1.0)
           );
    return  polya, polyb;
end


function  zig_zag_x( poly::Polygon2F, delta_x, delta_y, n::Int64,
                     shrink::Float64 = 0.95)

#    println( "n: ", n );
    for  i in 1:n
        p::Point2F = last( poly );
        q::Point2F = npoint( p[1] + delta_x, p[ 2 ] + delta_y );
        push!( poly, q  );
        delta_x = - shrink * delta_x;
    end
end


function  example_35()
    ( ! is_rebuild( "output/35" ) )  &&  return;

    delta = 0.01
    rho = 0.1
    h = 0.01;

    b = 2.0;
    
    P = Polygon2F() |> (-b, 0.0 ) |> (0.0, 0.0) |> (1.0, 0.0) |> (rho, h);

    zig_zag_x( P, 1.0 - 2.0*rho, h, 10, 1.0 );
    y = last( P )[ 2 ];
    P |> ( b, y + h );

    y = y + 2.0*h;
    

    Q = Polygon2F() |> (-b, y )  |> (0.0, y) |> (rho, y);
    zig_zag_x( Q, 1.0 - 2.0*rho, h, 10, 1.0 );
    y = last( Q )[ 2 ];
    Q |> (0.0, y ) |> (1.0, y + h );
    Q |> ( b, y + h );
    return  P, Q;
end



function  example_5()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.01, 6 );
    lp = last( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    push!( polyb,
           npoint( 0.0, 0.5 )
           ,
           npoint( max_x* ( 1-1.0/3.5 ), 0.5 )
           );
    zig_zag_x( polyb, -max_x/4.0, 0.02, 6 );
    pl = last( polyb );
    push!( polyb, npoint( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end


function  example_6_ext( zigs::Int64 )
    polya = Polygon2F( );
    polyb = Polygon2F( );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.001, 0.03, zigs );
    lp = last( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           npoint( 0.0, new_y ),
           npoint( max_x* 0.4, new_y ),
           npoint( max_x* 0.41, new_y + 0.5 ),
           npoint( max_x* 0.42, new_y ),
           npoint( max_x* ( 1-1.0/3.5 ), new_y )
           );
    zig_zag_x( polyb, -max_x/4.111111111, 0.025, zigs );
    pl = last( polyb );
    push!( polyb, npoint( 1.0*max_x, pl[ 2 ] ) );
    push!( polyb, npoint( 0.8*max_x, 1.1*pl[ 2 ] ) );
    push!( polyb, npoint( 1.0*max_x, 1.2*pl[ 2 ] ) );

    return  polya, polyb;
end

function  example_6()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.03, 16 );
    lp = last( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           npoint( 0.0, new_y ),
           npoint( max_x* 0.4, new_y ),
           npoint( max_x* 0.41, new_y + 0.5 ),
           npoint( max_x* 0.42, new_y ),
           npoint( max_x* ( 1-1.0/3.5 ), new_y )
           );
    zig_zag_x( polyb, -max_x/4.0, 0.025, 24 );
    pl = last( polyb );
    push!( polyb, npoint( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end



### Two segments crossing each other
function  example_7()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    xs = 4;

    push!( polya, npoint( 0.0, 0.1 ), npoint( xs*1.0, xs * 0.9 ) );
    push!( polyb,
           npoint( xs * 0.9, xs*1.0 )
           ,
           npoint( xs*0.7, 0 )
           );
    return  polya, polyb;
end


### Zig-zag on one curve, and straight segment on the other
function  example_8_ext( zigs::Int64 )
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/2.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.03, zigs );
    lp = last( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           npoint( 0.0, new_y ),
           npoint( max_x, new_y ) )

    return  polya, polyb;
end


function  example_9( zigs_a::Int64, zigs_b::Int64 )
    polya = Polygon2F( );
    polyb = Polygon2F( );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/1.1, 0 )
           );
    zig_zag_x( polya, -max_x/1.3, 0.03, zigs_a, 0.99 );
    lp = last( polya );
    pop!( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           npoint( 0.0, new_y ),
           npoint( max_x / 1.1, new_y )
           );
    zig_zag_x( polyb, -max_x/1.3, 0.025, zigs_b, 0.99 );
    pl = last( polyb );
    pop!( polyb );
    push!( polyb, npoint( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end


# Example 10
# Two zig-zag like curves that intersect in all their middle edges.
function  example_10( zigs_a::Int64, zigs_b::Int64 )
    P = Polygon2F( );
    Q = Polygon2F( );

    pref::Float64 = 0.3
    push!( P, npoint( pref/3.0, 0.0 ) )
    if  ( zigs_a > 1 )
        push!( P, npoint( pref, 0.0 ) );
    end
    push!( Q, npoint( 0.0, pref/3.0 ) );

    if  (zigs_a > 0 )  &&  ( zigs_a % 2 == 0 )
        zigs_a = zigs_a + 1
    end
    if  zigs_b % 2 == 0
        zigs_b = zigs_b + 1
    end
    d_a::Float64 = 0.4 / zigs_a;
    d_b::Float64 = 0.4 / zigs_b;

    f_down = false;
    y::Float64 = 0.0;
    for  i in 1:zigs_a
        f_down = !f_down;
        if  f_down
            y = 0.0
        else
            y = 1.0
        end
        push!( P, npoint( pref + d_a * i, y ) )
    end
    push!( P, npoint(1.0 - pref/3.0, 1.0 ) )

    f_left = false;
    x::Float64 = 0.0;
    for  i in 1:zigs_b
        f_left = !f_left;
        if  f_left
            x = 0.0
        else
            x = 1.0
        end
        push!( Q, npoint( x, pref + d_b * i ) )
    end
    push!( Q, npoint( 1.0, 1.0 - pref/3.0 ) )

    return  P, Q;
end

# Example 13: Based on example 7 - with additional fake vertices in
# the begining and end, to make the graph look nice...
function  example_13()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    xs::Float64 = 4.0;

    push!( polya,
           npoint( 0.0, 0.0 ),
           npoint( 0.1*xs, 0.0 ),
           npoint( 1.1*xs, xs * 1.0 ),
           npoint( 1.2*xs, xs * 1.0 ),
           );

    push!( polyb,
           npoint( 1.2 * xs, 0.0 ),
           npoint( 1.1 * xs, 0.0 ),
           npoint( 0.1 * xs, xs * 1.0 ),
           npoint( 0.0 * xs, xs * 1.0 )
           );
    return  polya, polyb;
end



function  example_14()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           npoint( 0.0, 0 )
           ,
           npoint( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.03, 16 );
    lp = last( polya );
    push!( polya, npoint( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           npoint( 0.0, new_y ),
           npoint( max_x, new_y ) );
#    pl = last( polyb );
 #   push!( polyb, npoint( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end

function   spiral_ext( P::Polygon2F, cen::Point2F, r::Float64, k::Int64,
                   shrink::Float64, range::Float64,
                   start_range::Float64 = 0.0 )::Polygon2F
    delta = range / k;

    for  i in 0:k
        alpha = start_range + delta * i;
        x = r * sin( alpha );
        y = r * cos( alpha );
        push!( P, npoint( cen[ 1 ] + x, cen[ 2 ] + y ) );
        r = r * shrink;
    end

    return  P;
end

function   spiral( cen::Point2F, r::Float64, k::Int64,
                   shrink::Float64, range::Float64 )::Polygon2F
    P = Polygon2F();
    spiral_ext( P, cen, r, k, shrink, range, pi / 2.0 );
    return  P;
end


function  example_15()
    max_x = 4.1;

    mid::Float64 = max_x/2.0;
    r = mid * 0.8;

    cen = npoint( mid, mid );
    polya = spiral( cen, 0.9*r, 6, 0.85, float(pi) );
    polyb = spiral( cen,     r,  15, 0.98, float(pi) );

    return  polya, polyb;
end


function  example_16()
    max_x = 4.1;

    mid::Float64 = max_x/2.0;
    r = mid * 0.8;

    cen = npoint( mid, mid );
    polya = spiral( cen, 0.9*r, 35*6, 0.99, 35*float(pi) );
    polyb = spiral( cen,     r,  40, 0.9, 3.999*float(pi) );

    return  polya, polyb;
end


function  example_17_dtw()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    k = 200;

    push!( polya, npoint(-1.0, -0.1 ), npoint( 6.0, -0.1 ) );

    spiral_ext( polya, npoint(2.0,0.3), 1.0, k, 1.0, Float64(pi),
                1.5 * Float64(pi) );

    push!( polyb, npoint(-1.0, 0.6 ) );
    spiral_ext( polyb, npoint(2.0,0.6), 1.1, k, 1.0, Float64(pi),
                1.5* Float64(pi) );

    return  polya, polyb;
end

### Two segments crossing each other
function  example_18()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    xs = 4;

    push!( polya,
        npoint( 0.0, 0.4999 * xs  ),
        npoint( 0.0, 0.499999 * xs  ),
        npoint( xs*1.0, xs * 0.5 ),
        npoint( xs*1.0, xs * 0.501 )
    );
    push!( polyb,
           npoint( xs * 0.4999, 0 ),
           npoint( xs * 0.49999, 0 )           ,
           npoint( xs * 0.600001, xs*1.0 ),
           npoint( xs * 0.6, xs*1.0 )
           );
    return  polya, polyb;
end


function  example_34()
    poly_a,poly_b = example_1();
    P = wiggle( poly_a, 40, 0.01 );
    Q = wiggle( poly_b, 40, 0.01 );
    return  P, Q;
end

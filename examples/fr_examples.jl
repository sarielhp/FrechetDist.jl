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

    push!( polya, point( 0.0, 0 ), point( xs*1.0, 0.0 ) );
    push!( polyb,
           point( 0.0, ys*0.1 ), point( xs * 0.45, ys * 0.1 ),
           point( xs*0.5, ys*0.9 ), point( xs*0.55  , ys*0.1),
           point( xs*1.0, ys*0.1 ) );
    return  polya, polyb;
end

function  example_2()
    polya = Polygon2F( );
    polyb = Polygon2F();

    push!( polya, point( 0.0, 0 ), point( 2.0, 0.0 ) );
    push!( polyb, point( 0.0, 1.0),  point( 0.4, 0.2 ),
           point( 0.6, 0.8 ),
           point( 1.0,0.3 ), point( 2.0, 0.3)  );
    return  polya, polyb;
end

function  example_3()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    push!( polya, point( 0.0, 0 ),
           point( 1.1, 0.0 ),
           point( 1.1, 0.1 ),
           point( 1.0, 0.1 ),
           point( 1.0, 0.2 ),
           point( 1.2, 0.2 ),
           point( 1.2, 0.0 ),
           point( 2.0, 0.0 ) );
    push!( polyb,
           point( 0.0, 0.3),
           point( 0.4, 0.3 ),
           point( 0.4, 0.6 ),
           point( 0.3, 0.6 ),
           point( 0.3, 0.7 ),
           point( 0.5, 0.7 ),
           point( 0.5, 0.3 ),
           point( 2.0, 0.3)  );
    return  polya, polyb;
end


function  example_4()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    push!( polya,
           point( 0.0, 0 )
           ,
           point( 0.0, 1 )
           );
    push!( polyb,
           point( 0.1, 0.0 )
           ,
           point( 1.1, 1.0)
           );
    return  polya, polyb;
end


function  zig_zag_x( poly::Polygon2F, delta_x, delta_y, n::Int64,
                     shrink::Float64 = 0.95)

#    println( "n: ", n );
    for  i in 1:n
        p::Point2F = last( poly );
        q::Point2F = point( p[1] + delta_x, p[ 2 ] + delta_y );
        push!( poly, deepcopy( q ) );
        delta_x = - shrink * delta_x;
    end
end

function  example_5()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           point( 0.0, 0 )
           ,
           point( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.01, 6 );
    lp = last( polya );
    push!( polya, point( max_x, lp[ 2 ] ) );

    push!( polyb,
           point( 0.0, 0.5 )
           ,
           point( max_x* ( 1-1.0/3.5 ), 0.5 )
           );
    zig_zag_x( polyb, -max_x/4.0, 0.02, 6 );
    pl = last( polyb );
    push!( polyb, point( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end


function  example_6_ext( zigs::Int64 )
    polya = Polygon2F( );
    polyb = Polygon2F( );

    max_x = 4.1;

    push!( polya,
           point( 0.0, 0 )
           ,
           point( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.001, 0.03, zigs );
    lp = last( polya );
    push!( polya, point( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           point( 0.0, new_y ),
           point( max_x* 0.4, new_y ),
           point( max_x* 0.41, new_y + 0.5 ),
           point( max_x* 0.42, new_y ),
           point( max_x* ( 1-1.0/3.5 ), new_y )
           );
    zig_zag_x( polyb, -max_x/4.111111111, 0.025, zigs );
    pl = last( polyb );
    push!( polyb, point( 1.0*max_x, pl[ 2 ] ) );
    push!( polyb, point( 0.8*max_x, 1.1*pl[ 2 ] ) );
    push!( polyb, point( 1.0*max_x, 1.2*pl[ 2 ] ) );

    return  polya, polyb;
end

function  example_6()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           point( 0.0, 0 )
           ,
           point( max_x/3.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.03, 16 );
    lp = last( polya );
    push!( polya, point( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           point( 0.0, new_y ),
           point( max_x* 0.4, new_y ),
           point( max_x* 0.41, new_y + 0.5 ),
           point( max_x* 0.42, new_y ),
           point( max_x* ( 1-1.0/3.5 ), new_y )
           );
    zig_zag_x( polyb, -max_x/4.0, 0.025, 24 );
    pl = last( polyb );
    push!( polyb, point( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end



### Two segments crossing each other
function  example_7()
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    xs = 4;

    push!( polya, point( 0.0, 0.1 ), point( xs*1.0, xs * 0.9 ) );
    push!( polyb,
           point( xs * 0.9, xs*1.0 )
           ,
           point( xs*0.7, 0 )
           );
    return  polya, polyb;
end


### Zig-zag on one curve, and straight segment on the other
function  example_8_ext( zigs::Int64 )
    polya = Polygon2F(  );
    polyb = Polygon2F(  );

    max_x = 4.1;

    push!( polya,
           point( 0.0, 0 )
           ,
           point( max_x/2.5, 0 )
           );
    zig_zag_x( polya, -max_x/4.0, 0.03, zigs );
    lp = last( polya );
    push!( polya, point( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           point( 0.0, new_y ),
           point( max_x, new_y ) )

    return  polya, polyb;
end


function  example_9( zigs_a::Int64, zigs_b::Int64 )
    polya = Polygon2F( );
    polyb = Polygon2F( );

    max_x = 4.1;

    push!( polya,
           point( 0.0, 0 )
           ,
           point( max_x/1.1, 0 )
           );
    zig_zag_x( polya, -max_x/1.3, 0.03, zigs_a, 0.99 );
    lp = last( polya );
    pop!( polya );
    push!( polya, point( max_x, lp[ 2 ] ) );

    new_y = lp[ 2 ] + 0.3;
    push!( polyb,
           point( 0.0, new_y ),
           point( max_x / 1.1, new_y )
           );
    zig_zag_x( polyb, -max_x/1.3, 0.025, zigs_b, 0.99 );
    pl = last( polyb );
    pop!( polyb );
    push!( polyb, point( max_x, pl[ 2 ] ) );

    return  polya, polyb;
end


# Example 10
# Two zig-zag like curves that intersect in all their middle edges.
function  example_10( zigs_a::Int64, zigs_b::Int64 )
    P = Polygon2F( );
    Q = Polygon2F( );

    pref::Float64 = 0.3
    push!( P, point( pref/3.0, 0.0 ), point( pref, 0.0 ) );
    push!( Q, point( 0.0, pref/3.0 ) );

    if  zigs_a % 2 == 0
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
        push!( P, point( pref + d_a * i, y ) )
    end
    push!( P, point(1.0 - pref/3.0, 1.0 ) )

    f_left = false;
    x::Float64 = 0.0;
    for  i in 1:zigs_b
        f_left = !f_left;
        if  f_left
            x = 0.0
        else
            x = 1.0
        end
        push!( Q, point( x, pref + d_b * i ) )
    end
    push!( Q, point( 1.0, 1.0 - pref/3.0 ) )

    return  P, Q;
end

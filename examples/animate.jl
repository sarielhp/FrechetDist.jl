#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist, FrechetDist.cg.polygon, Plots

# t ∈ [-1,1]
#heart(t) = 

#P = Polygon_fill( Polygon2F(), x->( sin(x), cos(x) )    , 0:0.01:2*pi )
#Q = Polygon_fill( Polygon2F(), x->( sin(2.0*x), cos(x) ), 0:0.01:2*pi )


# t ∈ [-1, 1]
heart(t) = ( 3.0*sin(t)* cos(t) * log( abs( t ) ),
             3.0*abs(t)^0.3*(cos(t))^0.5 -1.5 )


uheart(t) = ( 3.0*sin(t)* cos(t) * log( abs( t ) ),
             -3.0*abs(t)^0.3*(cos(t))^0.5 +1.5 )

heart(t) = (-16*(sin(t))^3,
            13.0*cos(t)-5.0*cos(2*t) -2.0 *cos(3*t) - cos(4.0*t) )
uheart(t) = (16*(sin(t))^3,
            -(13.0*cos(t)-5.0*cos(2*t) -2.0 *cos(3*t) - cos(4.0*t) ) )

P = Polygon_fill( Polygon2F(), heart, -pi:0.02:pi );
Q = Polygon_fill( Polygon2F(), uheart, 0:0.02:2.0*pi );

mrp = frechet_c_compute(P, Q);

M = Morphing_sample_uniformly( mrp, 200 );

anim = @animate for  i in 1:size(M,1)
    default( grid=:false, legend=false, framestyle=:none,
             label=:none, ticks=false, showaxis=false, lw=4 );
    plt = plot( 300, 300, 0, dpi = 200 );
    plot!( slice(P, 1), slice( P, 2 ),  lc=:red );
    plot!( slice(Q, 1), slice( Q, 2 ),  lc=:green );

    t = [ M[i,1:2]' ; M[i,3:4]' ];
    plot!( t[ :, 1 ], t[ :, 2 ], lw=8, lc=:blue );
end

gif(anim, "frechet.gif", fps = 10)

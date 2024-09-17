#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist, FrechetDist.cg Plots

P = Polygon_fill( Polygon2F(), x->( sin(x), cos(x) )    , 0:0.01:2*pi )
Q = Polygon_fill( Polygon2F(), x->( sin(2.0*x), cos(x) ), 0:0.01:2*pi )

mrp = frechet_c_compute(P, Q);

M = Morphing_sample_uniformly( mrp, 200 );

m_P = Polygon_as_matrix( P )
m_Q = Polygon_as_matrix( Q )

anim = @animate for  i in 1:size(M,1)
    plt = plot( 300, 300, 0,
                grid=:false, legend=false, framestyle=:none,
                label=:none, ticks=false, showaxis=false,
                dpi = 200 );
    plot!(plt, m_P[1,:], m_P[2,:], linewidth=4, legend=false,lc=:red );
    plot!(plt, m_Q[1,:], m_Q[2,:], lw=4, legend=false,lc=:green );
    t = [ M[i,1:2]' ; M[i,3:4]' ];
    plot!(plt, t[ :, 1 ], t[ :, 2 ], legend=false,lw=8, lc=:blue );
end

gif(anim, "frechet.gif", fps = 10)

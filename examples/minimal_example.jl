# /bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon

include( "graphics.jl" )


function  min_example()
    P = Polygon2F();
    Q = Polygon2F();

    push!( P, npoint( 0.0, 1.0 ), npoint( 10.0, 1.0 ), npoint( 10.0, 2.0 ),
           npoint( 20.0, 2.0 ) );
    push!( Q, npoint( 0.0, -1.0 ), npoint( 10.0, -2.0 ), npoint( 20.0, 3.0 ) );

    return  frechet_c_compute( P, Q );
end


mw = min_example()
output_morphing( mw, "morphing.pdf" );

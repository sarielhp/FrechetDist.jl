#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using Luxor

P = Polygon_random( 2, Float64, 10 );
for i ∈ eachindex( P )
    println( P[ i ]  );
end

m = zeros( length(P), 4 );
println( "=========================" );
for  i ∈ eachindex(P)
    #println( p );

    p = P[ i ];
    #println( typeof( x  ) );
    m[ i, : ] = [ p...; p... ];
end

println( m );
println( typeof( m ));
println( size( m ));

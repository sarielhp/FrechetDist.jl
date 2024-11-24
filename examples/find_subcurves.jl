#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using Graphs

struct  SubFrechetMatcher{N,T}
    X::Polygon{N,T};
    Y::Polygon{N,T};
    Delta::Int64;
    nX::Int64  # |X|
    nY::Int64  # |Y|
    G::SimpleDiGraph{Int64}

    function  SubFrechetMatcher( _X::Polygon{N,T},
                         _Y::Polygon{N,T},
                         _Delta::Int64 ) where {N,T}
        nX = cardin( _X );
        nY = cardin( _Y );
        
        nv = getNVerts( nX, nY, _Delta );
        return  new{N,T}( _X, _Y, _Delta, nX, nY,
                          SimpleDiGraph{Int64}( nv ) );
    end
end


function  SFM_compute_costs( M::SubFrechetMatcher{N,T} ) where {N,T}
    nX = M.nX;
    nY = M.nY;
    costs = Vector{Float64}();
    for  x in 1:nX
        for y in 1:nY
            push!( costs, Dist( M.X[ x ], M.Y[ y ] ) );
        end;
    end;
    sort!( costs );
    unique!( costs );
    return  costs;
end

function  dM(M::SubFrechetMatcher{N,T}, x::Int64, y::Int64 ) where  {N,T}
    return Dist( M.X[ x ], M.Y[ y ] )
end

function  SFM_build_graph( M::SubFrechetMatcher{N,T},
                         t::T ) where {N,T}
    nv = getNVerts( M );
    G = SimpleDiGraph{Int64}( nv );

    nX = M.nX;
    nY = M.nY;
    Delta = M.Delta;
    for  x::Int64 in 1:nX
        for y::Int64 in 1:nY
            # (x,y,1) ==> (x+1, y,1)
            if  (x + 1) <= nX
                c = max( dM( M, x, y ), dM( M, x+1, y ) );
                if  c <= t
                    add_edge!( G, M:(x,y,1), M:(x+1,y,1) );
                end
            end
            # (x,y,1) ==> (x+1, y,1)
            if  (y + 1) <= nY
                c = max( dM( M, x, y ), dM( M, x, y+1 ) );
                if  c <= t
                    add_edge!( G, M:(x,y,1), M:(x,y+2,1) );
                end
            end

            # (x,y,2) ==> (x, y,1)
            add_edge!( G, M:(x,y,2), M:(x,y,1) );

            # (x,y,1) ==> (x-Delta, y,2)
            add_edge!( G, M:(x,y,2), M:(x - Delta,y,2) );
        end;
    end;

    low = 1-Delta ;
    for  x::Int64 in low:nX
        for y::Int64 in 1:nY
            # (x,y,2) ==> (x-1, y,2)
            if  (x-1) >= low
                add_edge!( G, M:(x,y,2), M:(x-1,y,2) );
            end
            # (x,y,2) ==> (x, y-1,2)
            if  (y-1) >= 1
                add_edge!( G, M:(x,y,2), M:(x,y-1,2) );
            end
            # (x,y,2) ==> (x+1, y-1, 2)
            if  ( ( ( y-1 ) >= 1 )  &&  ( x+1 <= nX ) )
                add_edge!( G, M:(x,y,2), M:(x+1,y-1,2) );
            end
        end;
    end;

    return  G;
end

function  getNVerts( nX, nY, Delta )
    return  nX * nY  +  nY * ( nX + Delta );
end

function  getNVerts( M::SubFrechetMatcher{N,T} ) where {N,T}
    return  getNVerts( M.nX, M.nY, M.Delta );
end




function  get_vert( M::SubFrechetMatcher{N,T}, i::Int64 ) where {N,T}
    x::Int64 = 0;
    y::Int64 = 0;
    
    if  ( i <= M.nY * M.n )
        y = 1 + floor(Int64, (i-1)/ M.n)
        x = i - (y-1) * M.n;
        return  (x, y, 1 );
    end

    i -= M.nY * M.nX;
    len = M.nX + M.Delta;
    y = 1 + floor(Int64, (i-1) / len)
    x = i - (y-1) * len;
    return  ( x, y, 2 );
end

function (::Colon)(M::SubFrechetMatcher{N,T},
                    i::Tuple{Int64, Int64, Int64} ) where {N,T}
    t = i[3];
    x = i[1]; ## 1..n
    y = i[2]; ## 1..m
    if  t == 1
        return  (y-1) * M.nX + x;
    end
    @assert( t == 2 );
    @assert( x <= M.nX );
    @assert( y <= M.nY );
    @assert( x >= -M.Delta + 1 );
    x = x + M.Delta;
    return  M.nX * M.nY + (y-1) * (M.nX + M.Delta) + x;
end


######################################################################

function frechet_comp( P::Polygon{D,T}, Q::Polygon{D,T}
     )::Int64  where {D,T}
    ratio::Float64 = 5.0;

    println( "|P|:", cardin( P ), "  |Q| :", cardin( Q ) );

    for  i in 1:10
        println( "approx( #", cardin(P ), ", ", cardin(Q), ")   approx: ",
                 ratio );
        m = frechet_c_approx( P, Q, ratio );
        if  ( m.ratio == 1.0 )
            return  0;
        end
        lb = m.leash / m.ratio;

        #println( "m.leash: ", m.leash );
        ratio = (ratio - 1.0) / 6.0 + 1.0; # min( m.ratio, 1.01 );
        ratio = min( ratio, 1.1 );
#        println( "ratio: ", ratio );
        if  ( ratio <= 11.0 )
            println( "frechet_c_compute..." );
            m = frechet_c_compute( P, Q );
#            println( "... done" );
            return  0;
        end
#        if  i > 2
#            println( "RATIO  ", i, " : ", ratio );
#        end
    end
    println( "UNDECIDED" );

    return  0;
end


function  test_files( fl_a, fl_b )
    println( "---------------------------------------------" );
    println( fl_a );
    println( fl_b );
    println( "---------------------------------------------" );
    P = Polygon_read_file( fl_a );
    Q = Polygon_read_file( fl_b );

    M = SubFrechetMatcher( P, Q, 10 );
    costs = SFM_compute_costs( M );
    SFM_build_graph( M, costs[ length(costs) >> 1 ] );
    
    #frechet_comp( P, Q );
    println( "Done..." );
    #println( "Sign: ", sgn, "\n\n\n\n\n\n" );
end



num_args = length( ARGS );


if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 1 );
end

println( "compare_two_curves.jl [file1] [file2]" );

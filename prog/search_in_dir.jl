
#_media####
push!(LOAD_PATH, pwd())
push!(LOAD_PATH, pwd()*"/src/")


#using BenchmarkTools
using Parameters
using StaticArrays
using Distributions;
using Cairo
using LinearAlgebra
using Printf


using FrechetDist
using FrechetDist.cg


function search_dir(path,key)
    #println( readdir( path ) );
    return  filter(x->occursin(key,x), readdir(path))
end

###########################################################################
### Main starts here
###########################################################################
function find_closest_in_directory( path::String, threshold,
                                    si::Int64, sj::Int64 = 0 )
    files=search_dir( path, "plt" );

#    println( "bogi" );
    polygons = Vector{Polygon2F}();#  polygons;
    for  x in files
        filename = path * x;
        print( "Reading: ", filename );
        poly_a =  Polygon_read_plt_file( filename );
        push!( polygons,  poly_a );
        println( " #: ", cardin( poly_a ) );
    end
    println( "bogi2" );

    loc_i::Int64 = -1;
    loc_j::Int64 = -1;

    if  ( si > 0 )  &&  ( sj > 0 )  &&  (si != sj )
        println( "Special pair: ", si, " , ", sj );
        m = frechet_c_compute( polygons[ si ], polygons[ sj ] );
        exit( -1 );
    end

    min_dist = -1;
    len = length( files );
    widths = zeros(Float64, len )

    for  i in 1:len
        widths[ i ] = frechet_width_approx( polygons[ i ] );
    end
    
    dists_lb = zeros(Float64, len, len )
    dists_ub = zeros(Float64, len, len )

    for  i in 1:len
        dists_ub[ i, i ] = Inf
    end
    for  i in 1:( len - 1 )
        for  j  in  (i+1):len
            P = polygons[ i ];
            Q = polygons[ j ];
#            - widths[ i ] - widths[ j ];
            d = max( Dist( first( P ), first( Q ) ),
                     Dist( last( P ), last( Q ) ) );
            lb = d
            ub = d + widths[ i ] + widths[ j ];
            @assert( ub > 0.0 );
            dists_lb[ i, j ] = lb;
            dists_lb[ j, i ] = lb;
            dists_ub[ i, j ] = ub;
            dists_ub[ j, i ] = ub;
        end
    end

    pos = findmin( dists_ub )[ 2 ]
    up_val = dists_ub[ pos ];
    println( "pos :" , pos );
    println( "up_val: ", up_val );
    for  i in 1:( len - 1 )
        for  j  in  (i+1):len
            if  dists_lb[ i, j ] > up_val
                continue
            end
            print( "i; ", i, " j: ", j, " : " );
            m = frechet_c_approx( polygons[ i ], polygons[ j ], 2.0 );
            dist = m.leash;
            println( "  : ", dist, "        [",
                     dists_lb[ i,j ], " ... ",
                     dists_ub[ i,j ], "]" );

            if  ( dist < up_val )
                up_val = dist;
            end
            if  ( min_dist < 0 )  ||  ( dist < min_dist )
                loc_i = i;
                loc_j = j;
                min_dist = dist;
            end
#            println( "Curr   ", loc_i, ", ", loc_j, "  d: ",  min_dist );
        end
    end

    if  ( loc_i >= 0 )
        println( "Best pair found..." );
        println( path*files[ loc_i ] );
        println( path*files[ loc_j ] );
    end

    exit( -1 );
end




num_args = length( ARGS );
if   num_args == 1
    path = ARGS[ 1 ];
    println( path );
    find_closest_in_directory( path, 0.01, 0, 0 );

    exit( 0 );
end


println( "Wrong number of parameters!" );

exit( -1 );


#################################################################
#path = "data/012/trajectory/";
#base_a, base_b = "20080927005735.plt", "20081106012359.plt";
#base_a, base_b = "20080927005735.plt", "20081023012638.plt";
################################################################
# Excellent example:
#    filename_a = "data/026/trajectory/20090118104122.plt"
#    filename_b = "data/026/trajectory/20090225104803.plt"
################################################################



###############################################################################
###############################################################################
###############################################################################

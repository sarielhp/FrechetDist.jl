push!(LOAD_PATH, pwd()*"/src/")

#using BenchmarkTools
using Parameters
using StaticArrays
using Distributions;
using LinearAlgebra

using FrechetDist
using FrechetDist.cg
#using cg

function test_files( fl_a, fl_b )
    poly_a = Polygon_read_plt_file( fl_a );
    poly_b = Polygon_read_plt_file( fl_b );

    println( "#P: ", cardin( poly_a ) );
    println( "#Q: ", cardin( poly_b ) );
    


    println( "Approx 400..." );
    @time m_aprx_400 = frechet_c_approx( poly_a, poly_b, 10.0 );
    @time m_aprx_400 = frechet_c_approx( poly_a, poly_b, 10.0 );
    println( "Ratio: ", m_aprx_400.ratio );
    
    println( "Approx 40..." );
    m_aprx_40 = frechet_c_approx( poly_a, poly_b, 40.0 );
    @time m_aprx_40 = frechet_c_approx( poly_a, poly_b, 40.0 );
    println( "Ratio: ", m_aprx_40.ratio );
    
    println( "Approx 2.0..." );
    @time m_aprx_2 = frechet_c_approx( poly_a, poly_b, 2.0 );
    println( "2.0: ", m_aprx_2.ratio );

    println( "\nApprox 1.1..." );
    @time begin
        m_aprx_1 = frechet_c_approx( poly_a, poly_b, 1.1 );
    end
    println( "1.1: ", m_aprx_1.ratio );
    
    println( "Exact..." );
    @time m_exact = frechet_c_compute( poly_a, poly_b );

    println( "ve_r..." );
    @time m_ve_r = frechet_ve_r_compute( poly_a, poly_b );

end



num_args = length( ARGS );
if   num_args == 2
    test_files( ARGS[1], ARGS[2] );
    exit( 0 );
end

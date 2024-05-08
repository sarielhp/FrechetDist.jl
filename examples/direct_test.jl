push!(LOAD_PATH, pwd()*"/src/")

#using BenchmarkTools
#using Parameters
#using StaticArrays
#using Distributions;
#using LinearAlgebra
using TimerOutputs

using FrechetDist
using FrechetDist.cg
using CSV, DataFrames
#using cg

CSV_FILENAME = "output/results.csv";

CL_INPUT_P = "P #";
CL_INPUT_Q = "Q #";
CL_DESC = "Description";
CL_APRX_1_001 = "Aprx 1.001";
CL_APRX_1_01 = "Aprx 1.01";
CL_APRX_1_1 = "Aprx 1.1";
CL_APRX_2 = "Aprx 2";
CL_APRX_4 = "Aprx 4";
CL_EXACT = "Exact";
CL_VE_RETRACT = "VE Retractable";


function test_files( df::DataFrame, fl_a, fl_b, desc = "" )
    local  tmo;
    function  snano( str )::String
        return  string( TimerOutputs.time(tmo[ str ] ) / 1.0e9 );
    end

    tmo = TimerOutput()

    poly_a = Polygon_read_plt_file( fl_a );
    poly_b = Polygon_read_plt_file( fl_b );

    println( "#P: ", cardin( poly_a ) );
    println( "#Q: ", cardin( poly_b ) );

    ub_iters = cardin( poly_a ) * cardin( poly_b ) * 2;
    
    
    println( "Forcing compilation..." );
    m_exact_tmp = frechet_c_compute( poly_a, poly_b );
    #println( "\n\n\n\n\n\n\n\n\n\n\n\n\n\n" );

    local m_aprx_4;

    
    println( "Approx 4..." );
    @timeit tmo CL_APRX_4 begin
        m_aprx_4 = frechet_c_approx( poly_a, poly_b, 10.0 );
    end
    #println( "Ratio: ", m_aprx_4.ratio );

    println( "Approx 2.0..." );
    @timeit tmo CL_APRX_2 begin
        m_aprx_2 = frechet_c_approx( poly_a, poly_b, 2.0 );
    end
    #println( "2.0: ", m_aprx_2.ratio );

    println( "Approx 1.1..." );
    @timeit tmo CL_APRX_1_1 begin
        m_aprx_1 = frechet_c_approx( poly_a, poly_b, 1.1 );
    end
    #println( "1.1: ", m_aprx_1.ratio );

    #
    println( "\nApprox 1.001..." );
    @timeit tmo CL_APRX_1_001 begin
        m_aprx_1_001 = frechet_c_approx( poly_a, poly_b, 1.01 );
    end
    println( "1.001: ", m_aprx_1_001.ratio );

    println( "\nApprox 1.01..." );
    @timeit tmo CL_APRX_1_01 begin
        m_aprx_1_01 = frechet_c_approx( poly_a, poly_b, 1.01 );
    end
    println( "1.01: ", m_aprx_1_01.ratio );


    println( "Exact..." );
    @timeit tmo CL_EXACT begin
        m_exact = frechet_c_compute( poly_a, poly_b );
    end

    println( "ve_r..." );
    f_do_ve_r::Bool = ( ub_iters < 1000000000 );
    if  ( f_do_ve_r )
        @timeit tmo CL_VE_RETRACT begin
            m_ve_r = frechet_ve_r_compute( poly_a, poly_b );
        end
    end

    #show(tmo)


    #println( "\n" );
    #println( "Exact: ", nano( "Exact" ) );
    #println( TimerOutputs.time(tmo[ "Exact" ] ) );
    #println( TimerOutputs.time(tmo[ CL_APRX_4 ] ) );

    # Add the data to the data frame...
    push!( df, fill( "", ncol( df ) ) );
    r = nrow( df );

    df[ r, CL_INPUT_P    ] = string( cardin( poly_a ) );
    df[ r, CL_INPUT_Q    ] = string( cardin( poly_b ) );
    df[ r, CL_DESC       ] = desc;
    df[ r, CL_APRX_1_001 ] = snano( CL_APRX_1_001 );
    df[ r, CL_APRX_1_01  ] = snano( CL_APRX_1_01 );
    df[ r, CL_APRX_1_1   ] = snano( CL_APRX_1_1 );
    df[ r, CL_APRX_2     ] = snano( CL_APRX_2 );
    df[ r, CL_APRX_4     ] = snano( CL_APRX_4 );
    df[ r, CL_EXACT      ] = snano( CL_EXACT );
    if  ( f_do_ve_r )
        df[ r, CL_VE_RETRACT ] = snano( CL_VE_RETRACT );
    else
        df[ r, CL_VE_RETRACT ] = "TOO LARGE";
    end
    #=
CL_APRX_1_1 = "Aprx 1.1";
CL_APRX_2 = "Aprx 2";
CL_APRX_4 = "Aprx 4";
CL_EXACT = "Exact";
 = "VE Retractable";
=#
    print( df );
end


function  add_col( df::DataFrame, strs... )
    for  x in strs
        df[ :, x ] = String[];
    end
end



function  get_data_frame( filename::String )
    local df;
    if  ( isfile( filename ) )
        return   CSV.read( filename, DataFrame, types=String );
    end

    df = DataFrame();
    
    add_col( df, CL_INPUT_P, CL_INPUT_Q,  CL_APRX_1_001,
        CL_APRX_1_01, CL_APRX_1_1, CL_APRX_2, CL_APRX_4, CL_EXACT,
        CL_VE_RETRACT, CL_DESC );

    return  df;
end

df = get_data_frame( CSV_FILENAME )

println( df );

num_args = length( ARGS );
if   num_args == 2
    test_files( df, ARGS[1], ARGS[2] );
    CSV.write( CSV_FILENAME, df );
    exit( -1 );
end
if   num_args == 3
    test_files( df, ARGS[1], ARGS[2], ARGS[ 3 ] );
    CSV.write( CSV_FILENAME, df );
    exit( -1 );
end

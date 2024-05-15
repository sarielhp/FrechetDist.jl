#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using FrechetDist.cg
using PrettyTables
#using cg

include( "table_latex.jl" );

CSV_FILENAME = "output/julia_results.csv";

CL_INDEX      = "Input";
CL_INPUT_P    = "P #";
CL_INPUT_Q    = "Q #";
CL_DESC       = "Description";
CL_APRX_1_001 = "≈1.001";
CL_APRX_1_01  = "≈1.01";
CL_APRX_1_1   = "≈1.1";
CL_APRX_2     = "≈2";
CL_APRX_4     = "≈4";
CL_EXACT      = "Exact";
CL_VE_RETRACT = "VER";

function  str_int_w_commas( n::Int64 )
    if  ( abs( n ) < 1000 )
        return  string( n )
    end
    s = str_int_w_commas( floor( Int64, n/1000 ) );
    return s * @sprintf( ",%03d",  n % 1000 )
end


function test_files( df::DataFrame, fl_a, fl_b, desc = "" )
    local  tmo;
    function  snano( str )::String
        seconds::Float64 = TimerOutputs.time(tmo[ str ] ) / 1.0e9;
        s = @sprintf( "%.3f", seconds );
        return  s; #string(  );
    end

    tmo = TimerOutput()

    poly_a = Polygon_read_plt_file( fl_a );
    poly_b = Polygon_read_plt_file( fl_b );

    println( "#P: ", str_int_w_commas( cardin( poly_a ) ) );
    println( "#Q: ", str_int_w_commas( cardin( poly_b ) ) );

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

#    println( "Approx 2.0..." );
#    @timeit tmo CL_APRX_2 begin
#        m_aprx_2 = frechet_c_approx( poly_a, poly_b, 2.0 );
#    end
    #println( "2.0: ", m_aprx_2.ratio );

    println( "Approx 1.1..." );
    @timeit tmo CL_APRX_1_1 begin
        m_aprx_1 = frechet_c_approx( poly_a, poly_b, 1.1 );
    end
    #println( "1.1: ", m_aprx_1.ratio );

    #
    println( "\nApprox 1.001..." );
    @timeit tmo CL_APRX_1_001 begin
        m_aprx_1_001 = frechet_c_approx( poly_a, poly_b, 1.001 );
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

    df[ r, CL_INDEX      ] = string( r );
    df[ r, CL_INPUT_P    ] = str_int_w_commas( cardin( poly_a ) );
    df[ r, CL_INPUT_Q    ] = str_int_w_commas( cardin( poly_b ) );
    df[ r, CL_DESC       ] = desc;
    df[ r, CL_APRX_1_001 ] = snano( CL_APRX_1_001 );
    df[ r, CL_APRX_1_01  ] = snano( CL_APRX_1_01 );
    df[ r, CL_APRX_1_1   ] = snano( CL_APRX_1_1 );
#    df[ r, CL_APRX_2     ] = snano( CL_APRX_2 );
    df[ r, CL_APRX_4     ] = snano( CL_APRX_4 );
    df[ r, CL_EXACT      ] = snano( CL_EXACT );
    if  ( f_do_ve_r )
        df[ r, CL_VE_RETRACT ] = snano( CL_VE_RETRACT );
    else
        df[ r, CL_VE_RETRACT ] = "---";
    end
    #=
CL_APRX_1_1 = "≈ 1.1";
CL_APRX_2 = "≈ 2";
CL_APRX_4 = "≈ 4";
CL_EXACT = "Exact";
 = "VE Retractable";
=#
    println( df );

    ha = LatexHighlighter((d,i,j)->i % 2 != 0,
                          ["cellcolor{lightgray}","texttt"])
    hb = LatexHighlighter((d,i,j)->i % 2 == 0,
                          ["texttt"])

    dfn = deepcopy( df );
    select!(dfn, Not( CL_DESC ) )
    
    iox = open("output/julia_results.tex", "w");
    pretty_table( iox,dfn, header = names( dfn ), backend = Val(:latex),
                  highlighters = (ha, hb));

    #show( iox, "text/latex", df );
    close( iox );

    #################################################33
    
    dfn = deepcopy( df );
    iox = open("output/julia_results_2.tex", "w");
    select!(dfn, Cols( CL_INDEX, CL_DESC ) )

    pretty_table( iox,dfn, header = names( dfn ), backend = Val(:latex),
                  highlighters = (ha, hb));

    #show( iox, "text/latex", df );
    close( iox );
end




df = get_output_table( CSV_FILENAME )

println( df );
println( "\n" );

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

# Basic code to create the latex table...


push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
#using FrechetDist
#using FrechetDist.cg
using PrettyTables
#using cg

CSV_FILENAME = "output/results.csv";

CL_INDEX      = "Input";
CL_INPUT_P    = "P #";
CL_INPUT_Q    = "Q #";
CL_DESC       = "Description";
CL_APRX_1_001 = "≈1.001";
#CL_APRX_1_01  = "≈1.01";
CL_APRX_1_1   = "≈1.1";
CL_APRX_2     = "≈2";
CL_APRX_4     = "≈4";
CL_EXACT      = "Exact";
CL_VE_RETRACT = "VER";



function  add_col( df::DataFrame, strs... )
    for  x in strs
        df[ :, x ] = String[];
    end
end


function  get_output_table( filename::String )
    local df;
    if  ( isfile( filename ) )
        return   CSV.read( filename, DataFrame, types=String );
    end

    df = DataFrame();

    add_col( df, CL_INDEX, CL_INPUT_P, CL_INPUT_Q,  CL_APRX_1_001,
             #CL_APRX_1_01,
             CL_APRX_1_1,
             #CL_APRX_2,
        CL_APRX_4, CL_EXACT,
        CL_VE_RETRACT, CL_DESC );

    return  df;
end

#! julia

push!(LOAD_PATH, pwd()*"/src/")

using TimerOutputs
using Printf
using CSV, DataFrames
using FrechetDist
using BasicCompGeometry
import BasicCompGeometry: write_plt
using BasicCompGeometry.ReadWrite: read_geolife_plt
using PrettyTables


######################################################################

function plt_convert(filename)
    #    println( "FILENAME: [",filename,"]" );
    if ! isfile(filename)
        println("Error: No file ", filename)
        exit(-1)
    end
    s = readline(filename)
    if s != "Geolife trajectory"
        return
    end
    p = read_geolife_plt(filename)
    println("Rewriting: ", filename)
    write_plt(filename, p)
end

function plt_convert_in_file(filename)
    lines = readlines(filename)
    for i = 1:length(lines)
        plt_convert(lines[i])
    end
end


num_args = length(ARGS);

if num_args == 0
    println("plt2plt.jl [file1] [file2] ... ")
    exit(-1)
end

for i = 1:num_args
    println("--- ", ARGS[i])
    plt_convert_in_file(ARGS[i])
end

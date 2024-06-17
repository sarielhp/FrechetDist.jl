#! /bin/julia



############################################################################
# We get the queries test data from the frechet_distance project...
############################################################################

function  frechet_dist_project_get( dir_raw )
    fd_dir = dir_raw * "frechet/";
    #println( "FDDIr : ", fd_dir );
    if  ( ! isdir( fd_dir ) )
        mkdir( fd_dir )

        url = "https://gitlab.mpi-klsb.mpg.de/anusser/frechet_distance";
        x = run( ignorestatus( `git clone $url $fd_dir` ) )
    end

    return  fd_dir;
end

function   frechet_dist_project_get_query_tests( fd_dir )
    copied::Int64 = 0;
    queries_dir = fd_dir * "test_data/decider_benchmark_queries/";
    for (root, dirs, files) in walkdir( queries_dir )
        for file in files
            if  ! occursin( ".txt", file )
                continue;
            end
            dest_dir = "";
            if  occursin( "sigspatial", file )
                dest_dir = "data/queries/sigspatial/";
            elseif  occursin( "geolife", file )
                dest_dir = "data/queries/geolife/";
            elseif  occursin( "characters", file )
                dest_dir = "data/queries/characters/";
            else
                println( "ERROR: ", file );
                exit( -1 );
            end
            mkpath( dest_dir );
            new_file = dest_dir * file;
            old_file = joinpath(root, file);
            if  isfile( new_file )
                continue;
            end
            cp( old_file, new_file );
            copied = copied + 1;
            #println( "New: ", new_file );
            #println( "OF: ", old_file );
        end
    end
    (copied > 0)  &&  println( copied, " files copied!" );

    println( "Query test files in:\n"*
             "\t data/queries/sigspatial/ \n" *
             "\t data/queries/geolife/ \n" *
             "\t data/queries/characters/ \n\n");

end

function  test_command( cmds... )
    for  cmd  in cmds
#        println( "cmd: ", cmd );
        if  ( Sys.which( cmd ) == nothing )
            println( "Error:\n\t"*
                "Required system commad [", cmd,"] not found.\n\n"*
                 "Please install, and try again!\n\n" );
            exit( -1 );
        end
    end
end

##############################################33

test_command( "bash", "python3", "wget", "octave", "unzip", "tar" );

dir_raw = "raw/";

if  (! isdir( dir_raw ) )
    mkpath( dir_raw );
end

fd_dir = frechet_dist_project_get( dir_raw )
println( fd_dir );

dir_download = dir_raw * "frechet/test_data/benchmark/";

fetch_script = "fetch_and_convert_data.sh";


#cd( dir_download, cmd )

dir_data = "data/";
dir_geolife = "Geolife Trajectories 1.3/";
dir_geolife_trg = "data/geolife/";
dir_chars_trg  = "data/characters/";

if  !isdir( dir_download * dir_geolife )
    cmd = dir_download  * fetch_script;
    cp( "examples/replacement/"*fetch_script, dir_download * fetch_script,
        force=true );
#    println( "!!!" );
    run(Cmd( `./$cmd` ) );
#    println( "!!!" );
end

if  ! isdir( dir_geolife_trg )
    mkpath( dir_data );
    mv( dir_download * dir_geolife * "data/", dir_geolife_trg )
end

if  ! isdir( dir_chars_trg )
    mv( dir_download * "characters/data/", dir_chars_trg );
end

dir_sigspatial = "data/sigspatial";
if  ! isdir( dir_sigspatial )
    mv( dir_download * "sigspatial/", dir_sigspatial );
end

frechet_dist_project_get_query_tests( fd_dir );

println( "\n\n" *
         "Downloaded/convertex/etc data into data/ subdirectory."
    );

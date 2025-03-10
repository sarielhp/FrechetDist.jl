#! julia

push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"../src/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon
using CSV, DataFrames
using Downloads
using Dates;


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


function     mkpath_no_err( dr )
    if  isdir( dr )
        return
    end
    mkpath( dr );
end

##############################################33
# Splice data file...


#COL_ID = 14;
COL_TIME = 3;


function  diff_in_seconds( x::DateTime, y::DateTime )::Int64
    dt_diff = Dates.CompoundPeriod( x - y );
    if  ( isempty( dt_diff.periods ) )
        return  0;
    end
    diff = convert( Second, dt_diff );
    return  diff.value;
end

function  get_id_end( df, start::Int64, col::Int64 )
#    println( "\n\n\n\n" );
    #println( "get_id_end" );
    rows = nrow( df );
    #println( rows, " ", start );
    if  ( rows < start )
        return  -1000;
    end
    id = df[ start, col ];
    pos::Int64 = start;
    dt_prev = DateTime( df[ start, COL_TIME ], dateformat"y-m-d H:M:S.s" );
 #   println( dt_prev );
#    println( typeof( dt_prev ) );
    seconds_in_day = 24 * 3600;
    gap = seconds_in_day / 2;  # was 3 * seconds_in_day
    while  pos <= rows
        if  ( df[ pos, col ] != id )
            return pos - 1;
        end
        dt_curr = DateTime( df[ pos, COL_TIME ], dateformat"y-m-d H:M:S.s" );
        diff::Int64 = diff_in_seconds( dt_curr, dt_prev );

        if  ( diff > gap )
 #           println( "Diff: ", diff, "    ",  diff/seconds_in_day );
            return  pos - 1;
        end

        dt_prev = dt_curr;
        pos = pos + 1;
    end
    return  pos - 1;
end

function  print_headers( df )
    nmx = names( df );
    cols = length( nmx );
    for  i in 1:cols
        println( "[", i, "] ", nmx[ i ], ": ",  df[ 1, i ] );
    end
end

function  splice( df, base_dir, id_col )
    pos::Int64 = 1;
    prev_name::String = "nope nope";
    count::Int64 = 1;
    #println( "pos: ", pos );
    while  pos > 0
        #println( "pos: ", pos );
        fin::Int64 = get_id_end( df, pos, id_col );
        if   fin < 0
            break;
        end
        poly_id = df[ pos, id_col ];
        poly = Polygon2F();
        for  i in  pos:fin
            x::Float64 = df[ i, 4 ];
            y::Float64 = df[ i, 5 ];
            #y = parse( Float64, df[ i, 4 ] )
            push!( poly, npoint( x, y ) );
        end
        poly_name = string( poly_id );
        if  ( prev_name != poly_name )
            count = 1;
        else
            count = count + 1;
        end
        #println( typeof( poly_name ) );
        out_name = string(poly_name) * "_" * string( count ) * ".plt";

        filename = base_dir * "/" * out_name;
        println( "Writing: ", filename );
        write_to_file( poly, filename );
        prev_name = poly_name;
        #if  ( fin > 0 )
        #println( pos, "...", fin, "  ", df[ pos, id_col ] );
        #    end
        pos = fin + 1;
    end
end

function print_columns( df )
    cols = ncol( df );
    for  i  in 1:cols
        println( i, ": ", names( df )[ i ] );
    end
end

function  splice_data_file( data_filename, dir_out )
    df = CSV.read( data_filename, DataFrame );

    id_col = findfirst(item -> item == "tag-local-identifier", names( df ) )


    splice( df, dir_out, id_col );
end



##############################################33
# Birds...


function   get_birds( dir_raw, dir_birds )
    if  isdir( dir_birds )
        return
    end

    dir_raw_birds = dir_raw * "birds/";
    mkpath_no_err( dir_raw_birds );

    birds_base = "eastern_flyway_spring_migration_of_adult"*
                 "_white_storks"*
                 "__data_from_rotics_et_al__2018__gps.csv";
    birds_name = dir_raw_birds * birds_base;
    birds_url = "https://datarepository.movebank.org/bitstreams"*
                "/ed3529e6-ce20-4237-8875-9d35fdbf9a0f/download";
    birds_alt_url = "https://frechet.xyz/download/";
    if  ! isfile( birds_name )
        println( "Downloading birds data... Might take a while..." );
        run( `wget "$birds_url"` )
        sz = filesize( "download" );
        println( "SIZE: ", sz );
        if  sz < 273153300
            rm( "download" );
            url = birds_alt_url * birds_base;
            run( `wget "$url"` )
            mv( birds_base, birds_name );
        else
            mv( "download", birds_name );
        end
        println( "\tDownloading done!" );
    end
    mkpath_no_err( dir_birds );

    splice_data_file( birds_name, dir_birds );
end


################################################################333
# Pigeons

function  get_pigeons_data( dir_raw, dir_pigeons )
    dir_raw_pigeons = dir_raw * "/pigeons/";
    mkpath_no_err( dir_raw_pigeons );
    url_base = "https://frechet.xyz/download/";
    base_name = "right_hemisphere_advantage_in_route_fidelity_in"*
                "_homing_pigeons__data_from_pollonara_et_al__2017__csv.zip"
    zip_name = dir_raw_pigeons*base_name;
    zip_csv_name = "Right hemisphere advantage in route fidelity in "*
               "homing pigeons (data from Pollonara et al. 2017).csv";
    csv_name = dir_raw_pigeons * "pigeons.csv";
    if  ! isfile( csv_name )
        Downloads.download( url_base*base_name, zip_name );
        run( Cmd( `unzip $zip_name -x '__MACOSX/*'` ) );
        mv( zip_csv_name, csv_name );
    end

    if  ! isdir( dir_pigeons ) 
        mkpath_no_err( dir_pigeons );
        splice_data_file( csv_name, dir_pigeons );
    end
end

##############################################33

test_command( "bash", "python3", "wget", "octave", "unzip", "tar" );

dir_raw = "raw/";

mkpath_no_err( dir_raw );

fd_dir = frechet_dist_project_get( dir_raw )
println( fd_dir );

dir_download = dir_raw * "frechet/test_data/benchmark/";

fetch_script = "fetch_and_convert_data.sh";
convert_script = "geolife_converter.py";


#cd( dir_download, cmd )
dir_raw = "raw/";
dir_data = "data/";
dir_geolife = "Geolife Trajectories 1.3/";
dir_geolife_trg = "data/geolife/";
dir_chars_trg  = "data/characters/";

if  !isdir( dir_download * dir_geolife )
    cmd = dir_download  * fetch_script;
    cp( "examples/replacement/"*fetch_script, dir_download * fetch_script,
        force=true );
    cp( "examples/replacement/"*convert_script, dir_download * convert_script,
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

println( "data/sigspatial..." );
dir_sigspatial = "data/sigspatial";
if  ! isdir( dir_sigspatial )
    mv( dir_download * "sigspatial/", dir_sigspatial );
end

println( "data/queries..." );
if  ! isdir( "data/queries" )
    frechet_dist_project_get_query_tests( fd_dir );
end

println( "data/birds..." );
dir_birds = "data/birds";
if  ! isdir( dir_birds )
   get_birds( dir_raw, dir_birds )
end

println( "data/pigeons..." );
dir_pigeons = "data/pigeons"
if  ! isdir( dir_pigeons )
   get_pigeons_data( dir_raw, dir_pigeons )
end

println( "\n\n" *
         "Downloaded/convertex/etc data into data/ subdirectory."
    );

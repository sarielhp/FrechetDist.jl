#! /bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using Downloads;
using FrechetDist
using FrechetDist.cg


tryusing(pkgsym) = try
    @eval using $pkgsym
    return true
catch e
    return e
end

function  print_progress(total::Integer, now::Integer)
    print( now, " / ", total, "      \r" );
end

#if tryusing(:UrlDownload) != true
#    using  Pkg
#    Pkg.add("UrlDownload")
#end

function  test_zip( filename )::Bool
    if  ( ! isfile( filename ) )
        return  false;
    end
    x = run( ignorestatus( `unzip -t $filename` ))
    if  x.exitcode == 0
        return true
    end

    rm( filename );

    return  false;
end

function  test_tgz( filename )::Bool
    if  ( ! isfile( filename ) )
        return  false;
    end
    x = run( ignorestatus( `tar -tf $filename` ))
    if  x.exitcode == 0
        return true
    end

    rm( filename );

    return  false;
end

######################################################################3

function  plt_convert( filename, dest )
    if  isfile( dest )
        return
    end
    if  ! isfile( filename )
        println( "Error: No file ", filename );
        exit( -1 );
    end
    s = readline( filename );
    if   s != "Geolife trajectory"
        return
    end
    p = Polygon_read_plt_orig_file( filename );
    print( " => ", dest, "           \r" );
    Polygon_write_to_file( p, dest );
end


function  geolife_get( dir_raw )
    geo_dir = dir_raw * "geolife_in/";
    if  ( isdir( dir_raw * "geolife_in/" ) )
        println( "Geolife in directory exists." *
            " Assuming download was successful" )
        return geo_dir;
    end
    geolife_zip = dir_raw * "geolife.zip";
    url_geolife = "https://download.microsoft.com/download/" *
                  "F/4/8/F4894AA5-FDBC-481E-9285-D5F8C4C4F039" *
                  "/Geolife%20Trajectories%201.3.zip";

    if  ( ! test_zip( geolife_zip ) )
        Downloads.download( url_geolife, geolife_zip,
            progress=print_progress  );
        if  ! test_zip( geolife_zip );
            println( "Error! Unable to download geolife data file." );
            exit( -1 );
        end
    end
    if  ( ! isdir( geo_dir ) )
        run( `unzip $geolife_zip` );
        mv( "Geolife Trajectories 1.3", dir_raw * "geolife_in/" );
    end

    return  geo_dir;
end

function  geolife_convert( geo_dir, data_dir )
    cnt::Int64 = 0;
    files_c::Int64 = 0;
    println( "" );
    for (root, dirs, files) in walkdir( geo_dir )
        for file in files
            if  ! occursin( ".plt", file )
                continue;
            end
            files_c = files_c + 1;
            rt = split.(root, "/")

            new_dir = data_dir * rt[ 4 ] * "/" *
                      lowercase( rt[5] ) * "/";
            new_file = new_dir * file;
            if  ( isfile( new_file ) )
                continue;
            end
            old_file = joinpath(root, file);

            if  ( ! isdir( new_dir ) )
                mkpath( new_dir );
            end
            cnt = cnt + 1;
            #println( old_file, "=>", new_file )
            plt_convert( old_file, new_file )
        end
    end
    if  ( cnt > 0 )
        println( "\n\n" );
    end

    println( data_dir, " created. ", cnt, "/", files_c, " converted." );
end


function  geolife_do( dir_raw )
    geo_dir = geolife_get( dir_raw );
    geolife_zip = dir_raw * "geolife.zip";
    geolife_convert( geo_dir, "data/geolife/" );

    println( "\nYou can safely delete: \n\t", geo_dir * " (directory)\n"*
         "\t", geolife_zip, "\n\n" );
end

############################################################################

function  sigspatial_get( dir_raw )

    sig_dir = dir_raw * "sigspatial/";
    if  ( isdir( sig_dir ) )
        println( "Geolife in directory exists." *
            " Assuming download was successful" )
        return sig_dir;
    end
    sigspatial_tgz = dir_raw * "sigspatial.tgz";
    url_sigspatial = "http://sfsp.mpi-inf.mpg.de/shortest-sf.tgz";

    if  ( ! test_tgz( sigspatial_tgz ) )
        Downloads.download( url_sigspatial, sigspatial_tgz,
            progress=print_progress  );
        if  ! test_tgz( sigspatial_tgz );
            println( "Error! Unable to download geolife data file." );
            exit( -1 );
        end
    end
    if  ( ! isdir( sig_dir ) )
        run( `tar xf $sigspatial_tgz` );
        mv( "files", sig_dir );
    end

    return  sig_dir;
end

function  sigspatial_do( dir_raw )
    println( "\n===== Sigspatial:\n" );
    if  isdir( "data/sigspatial" )
        println( "data/sigspatial/ already exists.\n" );
        return;
    end
    dir_sig = sigspatial_get( dir_raw );
    mv( dir_sig, "data/sigspatial/" );
end

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

##############################################33

dir_raw = "raw/";

if  (! isdir( dir_raw ) )
    mkpath( dir_raw );
end

geolife_do( dir_raw );

sigspatial_do( dir_raw );

fd_dir = frechet_dist_project_get( dir_raw )
frechet_dist_project_get_query_tests( fd_dir );

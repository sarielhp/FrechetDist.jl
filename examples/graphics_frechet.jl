# Graphics_frechet.jl

using Parameters
using Distributions;
using Cairo
using LinearAlgebra
using Printf
using Plots
using PrettyTables
using Dates

#----------------------------------------------------------------
# Output the morphing to a pdf file
function  output_morphing( m::Morphing{N,T}, filename )  where {N,T}
    c,cr,bb = cairo_setup( filename, [ m.P, m.Q ], true );

    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_line_width(cr, 1.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    P = m.P;
    Q = m.Q;
    PS, QS = Morphing_as_polygons( m );

    set_line_width(cr, 0.5);
    set_source_rgb(cr, 1.0, 0.0, 0.0 );

    nv::Int64 = cardin( PS );
    for  i in 1:nv
        #loc::Point2I  = sol.pnts[ i ];

        po::Point{N,T} = PS[ i ]
        qo::Point{N,T} = QS[ i ]

#        println( loc.x[ 1 ], ", ", loc.x[ 2] );

        move_to( cr,  po[1], po[2] )
        line_to( cr, qo[1], qo[2] );
        Cairo.stroke( cr );
    end

    set_line_width(cr, 1.0);
    set_source_rgb(cr, 0.0, 1.0, 0.0 );
    draw_polygon( cr, P );
    set_source_rgb(cr, 0.0, 0.0, 1.0 );
    draw_polygon( cr, Q );

    set_source_rgb( cr, 1,0,0 );

    Cairo.finish(c);
end


function  output_frechet_movie( m::Morphing{N,T},
                                filename::String,
    total_frames::Int64 = 800,
    f_show_vertices::Bool = false ) where {N,T}

    if  isfile( filename )
        println( "\n\n",  filename, " already exists...\n\n\n" );
        return;
    end

    c,cr,bb = cairo_setup( filename, [ m.P, m.Q ], true );

    #   Cairo.save( cr );

    set_line_width(cr, 0.5);
    set_source_rgb(cr, 1.0, 0.0, 0.0 );


    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    if   np != nq
        println( "Error np!= nq" );
        exit( -1 );
    end

    lpout::Vector{Float64} = Polygon_prefix_lengths( pout )
    lqout::Vector{Float64} = Polygon_prefix_lengths( qout )

    lens::Vector{Float64} = lpout + lqout;

    steps = length( lens ) - 1;
    total_length = last( lens );
    #println( lens );


    # compute how many frames for each leg...
    frames = compute_frames( pout, qout, total_frames )

    total_output_frames = sum( frames );
    println( "Total number of frames: ", total_output_frames );
    skip::Int64 = max( 1, floor( total_output_frames / total_frames ) );
    f_do_skip::Bool = true ;

    for  i in 1:steps
        count = 0;
        if  frames[ i ] == 0
            continue;
        end
        if  ( ! f_do_skip )
            f_output_frame = true;
        else
            f_output_frame = ( count % skip == 0 )  ||  ( i == steps );
        end
        count = count + 1;
        if  ( f_output_frame )
            sp = Segment( pout[ i ], pout[ i + 1 ] );
            sq = Segment( qout[ i ], qout[ i + 1 ] );

            draw_frames( cr, sp, sq, frames[ i ], m.P, m.Q, bb );
        end
    end

    Cairo.finish(c);
end

function  rmx( tmp_filename )
    if isfile( tmp_filename )
        rm( tmp_filename );
    end
end


function  output_frechet_movie_mp4( m::Morphing{N,T},
                                filename::String,
    total_frames::Int64 = 800,
    f_show_vertices::Bool = false
) where {N,T}
    f_debug::Bool = false;
    cm = ContextMovie(BBox2F(), BBox2F(), 0, 0, 0, "/tmp/r_draw/",
        f_show_vertices );
    cm.bb, cm.bbo = compute_bounding_boxes( [ m.P, m.Q ] );
    cm.iheight, cm.iwidth = get_image_dims( cm.bbo );

    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    if   np != nq
        println( "Error np!= nq" );
        exit( -1 );
    end

    #set_transform( cr, iwidth, bbo );

    #set_transform( cr, iwidth, bbo );

    #c = Cairo.CairoPDFSurface( filename, iwidth, iheight );
    #cr = CairoContext(c);
    # Create temporary directory for images...
    if  isdir( cm.dir )
        rm( cm.dir, recursive=true)
    end
    mkdir( cm.dir );


    frames, steps = compute_frames( pout, qout, total_frames )
    total_output_frames = sum( frames );
    skip::Int64 = max( 1, floor( total_output_frames / total_frames ) );
    f_do_skip::Bool = true ;

    vec_rf = Vector{RecFrame}();
    frame_count::Int64 = 0;

    # We first calculate the frames we need... into vec_rf...
    for  i in 1:steps
        count = 0;
        if  frames[ i ] == 0
            continue;
        end
        if  ( ! f_do_skip )
            f_output_frame = true;
        else
            f_output_frame = ( count % skip == 0 )  ||  ( i == steps );
        end
        count = count + 1;

        if  ( f_output_frame )
            sp = Segment( pout[ i ], pout[ i + 1 ] );
            sq = Segment( qout[ i ], qout[ i + 1 ] );

            draw_frames_images( cm, sp, sq, frames[ i ], m.P, m.Q, vec_rf );
        end
    end


    f_debug  &&  println( "Splitting to threads... Sit back ;)" );
    f_debug  &&  println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)

    ####
    # We now call ffmpeg to create the movie...
    # ffmpeg -r 10 -i temp/ do.mp4

#    cmd = ( " -r 10 -i " *   "
#          * filename );
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
#    println( "ffmpeg $options" );
    f_debug && println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               cm.dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);

    f_debug &&  println( "Rencoding with handbrake..." );

    # HandBrakeCLI -Z  -i movie.mp4  -o movie_2.mp4
    options_2 = [ "-Z", "Android 1080p30", "-i", tmp_filename,
                 "-o", filename ];
    output_2 = read(pipeline( `HandBrakeCLI $options_2`,
                            stderr="/tmp/errs_2.txt" ), String);
    rmx( tmp_filename );
    if  isfile( filename )
        println( "Created: ", filename );
    end
end




function  frames_generate_o( cm::ContextMovie, P, Q, vec_rf )
    for rf in vec_rf
        draw_image_frame_o( cm::ContextMovie, P, Q, rf )
    end
    return  1
end

function  output_ortho_frechet_movie_mp4( m::Morphing{N,T},
                                          filename::String,
                                          total_frames::Int64 = 200,
                                          ) where {N,T}
    cm = ContextMovie(BBox2F(), BBox2F(), 0, 0, 0, "/tmp/r_draw/",
                      false  );
    cm.bb, cm.bbo = compute_bounding_boxes( [ m.P, m.Q ] );
    cm.iheight, cm.iwidth = get_image_dims( cm.bbo );

    pout,qout = Morphing_as_polygons( m );

    np::Int64 = cardin( pout );
    nq::Int64 = cardin( qout );
    @assert( np == nq );

    if  isdir( cm.dir )
        rm( cm.dir, recursive=true)
    end
    mkdir( cm.dir );

    vec_rf = Vector{ORecFrame}();

    delta = 1.0 / ( total_frames - 1 );
    t::Float64 = 1.0;

    # We first calculate the frames we need... into vec_rf...
    for  i in 1:total_frames
        R = Polygon_convex_comb( pout, qout, t );

        ofr = ORecFrame( i, t, R );
        push!( vec_rf, ofr );
        t = t - delta;
    end

    #println( "Splitting to threads... Sit back ;)" );

    #println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate_o( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)

    encode_to_mp4( cm.dir, filename );
end


#P::Polygon{N,T}, Q::Polygon{N,T}, Pe, Qe,
function  output_frechet_diagram( m::Morphing{N,T}, filename )  where {N,T}

    P_coords::Vector{T} = get_diagram_locs( m.pes, m.P );
    Q_coords::Vector{T} = get_diagram_locs( m.qes, m.Q );

    len = length( P_coords );
    poly = Polygon2F();
    for  i in 1:len
        push_smart!( poly, npoint( P_coords[ i ], Q_coords[ i ] ) );
    end

    psum::Vector{Float64} = Polygon_prefix_lengths( m.P )
    qsum::Vector{Float64} = Polygon_prefix_lengths( m.Q )

    c,cr,bb = cairo_setup( filename, [ poly ], true );

    set_source_rgb(cr,0.9,0.0,0.0);
    set_line_width(cr, 1.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    plen = cardin( m.P );
    qlen = cardin( m.Q );
    ymax = last( qsum );
    xmax = last( psum );

    for  i in 1:plen
        move_to( cr, psum[i], 0 )
        line_to( cr, psum[i], ymax );
        Cairo.stroke( cr );
    end

    for  i in 1:qlen
        move_to( cr,    0, qsum[i] );
        line_to( cr, xmax, qsum[i] );
        Cairo.stroke( cr );
    end
    set_line_width(cr, 1.00);
    set_source_rgb(cr, 0.0, 1.0, 0.0 );
    draw_polygon( cr, poly );
    Cairo.finish(c);
end


function  get_diagram_locs( PE::Vector{EventPoint{N,T}}, P::Polygon{N,T}
                            ) where {N,T}

    prefixes::Vector{Float64} = Polygon_prefix_lengths( P )

    pout = Vector{T}()

    len = length( PE );
    push!( pout, 0 );

    i = 2;
    while  ( i <= len )
        ep = PE[ i ];
        if  ep.type == PT_VERTEX
            push!( pout, prefixes[ ep.i ] );
            i = i + 1;
            continue;
        end

        loc = ep.i;
        edge_length = prefixes[ loc + 1 ] - prefixes[ loc ];
        push!( pout, prefixes[ loc ] + ep.t * edge_length );

        i = i + 1
    end

    return  pout
end



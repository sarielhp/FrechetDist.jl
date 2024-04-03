# Originally contributed by S. Har-Peled
# under MIT License

push!(LOAD_PATH, pwd()*"/src/")


#using BenchmarkTools
using Parameters
using StaticArrays
using Distributions;
using Cairo
using LinearAlgebra
using Printf
using Plots
using PrettyTables
using Dates

using FrechetDist
using FrechetDist.cg
#using cg

include( "fr_examples.jl" )


FrechetStr::String = "Fréchet";


function  is_mkdir( dir )
    if  ! isdir( dir )
        mkdir( dir );
    end
end


VecFloat = Vector{Float64};
VecVecFloat = Vector{VecFloat};

function  draw_hippodrome( cr, p::Point2F, q::Point2F, r::Float64 )
    r = r * 4.0;

    arc( cr, p[1], p[2], r, 0.0, 2.0 * pi );
    fill_preserve( cr );
    arc( cr, q[1], q[2], r, 0.0, 2.0 * pi );
    fill_preserve( cr );

    v = (p - q) / Dist( p, q );
    u = point( -v[2], v[1] );
    px = p - r * u;
    qx = q - r * u;
    qy = q + r * u;
    py = p + r * u;

    move_to(cr, px[ 1 ], px[ 2 ] );
    line_to(cr, qx[ 1 ], qx[ 2 ] );
    line_to(cr, qy[ 1 ], qy[ 2 ] );
    line_to(cr, py[ 1 ], py[ 2 ] );
    line_to(cr, px[ 1 ], px[ 2 ] );
    close_path(cr);
    stroke_preserve(cr);
    Cairo.fill(cr);
    Cairo.stroke( cr );
end

function  draw_polygon_w_offs( cr, P::Polygon2F, offs::VecFloat )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        p = P.pnts[ i - 1 ];
        #p = po.x;
        q = P.pnts[ i ];
        #q = qo.x;

        r = max( offs[ i - 1 ], offs[ i ] );
        #println( " r: ", r );

        set_source_rgb(cr, 1.0, 0.0, 0.8);
        #arc( cr, p[1], p[2], r, 0.0, 2.0 * pi );
        #fill_preserve( cr );

        draw_hippodrome( cr, p, q, r );

        #Cairo.stroke( cr );
        #cr->set_source_rgba(0.0, 0.0, 0.8, 1.0);
        #move_to( cr,  p[1], p[2] )
        #line_to( cr, q[1], q[2] );
        #println( q[1], " ", q[2] );
        Cairo.stroke( cr );
    end
end


function  draw_polygon( cr, P )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        po = P.pnts[ i - 1 ];
        p = po.x;
        qo = P.pnts[ i ];
        q = qo.x;
        move_to( cr,  p[1], p[2] )
        line_to( cr, q[1], q[2] );
        #println( q[1], " ", q[2] );
    end
    Cairo.stroke( cr );
end

function  draw_polygon_vertices( cr, P, r::Float64 )
    nv::Int64 = cardin( P );
    #r = Polygon_length( P ) / (100.0*cardin( P ));
    for  i in 1:nv
        p = P.pnts[ i ];
#        set_line_width(cr, 2.00);
#        Cairo.set_source_rgb( cr, 0, 0, 0);
#        println( "RRR = ", r );
        Cairo.arc( cr, p[1], p[2], r, 0.0, 2 * pi);
        Cairo.fill(cr);
#        Cairo.arc( cr, p[1], p[2], r, 0.0, 2 * pi);
#        Cairo.set_source_rgb( cr, 1.0, 1.0, 0);
    end
    Cairo.stroke( cr );
end

function  draw_bbox( cr, bb, scale )
    pa = BBox_bottom_left( bb );
    pc = BBox_top_right( bb );

    pb = point( pc[1], pa[2] );
    pd = point( pa[1], pc[2] )

    pa = pa * scale;
    pb = pb * scale;
    pc = pc * scale;
    pd = pd * scale;
    ;
    move_to( cr, pa[1], pa[2] )
    line_to( cr, pb[1], pb[2] );
    line_to( cr, pc[1], pc[2] );
    line_to( cr, pd[1], pd[2] );
    line_to( cr, pa[1], pa[2] );

    Cairo.stroke( cr );
end


function  compute_bounding_boxes( list::VecPolygon2F )
    bb::BBox2F = BBox2F();

    BBox_bound( bb, list );
    BBox_expand( bb, 1.05 );
    bbo::BBox2F = deepcopy( bb );
    BBox_expand( bbo, 1.05 );

    return  bb, bbo
end

function  get_image_dims( bbo )

    width::Float64 = 1024.0;
    theight::Float64 = 0.0;

    while ( true )
        theight = width * BBox_width( bbo, 2 ) / BBox_width( bbo, 1 );
        if  theight < 2048.0
            break;
        end

        width = width / 2.0;
    end

    iheight::Int64 = convert( Int64, 16 * ceil( theight / 16 ) )
    iwidth::Int64 = convert( Int64, 16 * ceil( width / 16 ) )

    return  iheight,iwidth;
end

function  set_transform( cr, iwidth::Int64, iheight::Int64,
                         bbo::BBox2F )
    xcal = convert( Float64, iwidth) / BBox_width( bbo, 1 );

#        ycal = convert( Float64, iheight) / BBox_width( bbo, 2 );
#    if  ( ( (xcal / 5.0 ) < ycal ) && ( ycal <= (xcal * 5.0 ) ) )
#        ycal = xcal;
#    end

    Cairo.scale( cr, xcal, xcal );
    bl = BBox_bottom_left( bbo );
    Cairo.translate( cr, -bl[ 1 ], -bl[ 2 ]);
end

function  cairo_setup( filename::String, list::VecPolygon2F,
                       f_pdf::Bool = true )
    bb, bbo = compute_bounding_boxes( list );
    iheight, iwidth = get_image_dims( bbo );

#    set_transform( cr, iwidth, bbo );

    local c
    if (  f_pdf )
        c = Cairo.CairoPDFSurface( filename, iwidth, iheight );
    else
        c = CairoRGBSurface( iwidth, iheight );
    end
    cr = CairoContext(c);

    set_transform( cr, iwidth, iheight, bbo );

    if  ( ! f_pdf )
        set_source_rgb(cr, 1, 1, 1);
        paint(cr);
    end

    return  c,cr,bb;
end


function  output_polygons_to_file(  list::VecPolygon2F, filename,
                                    f_pdf::Bool,
                                    f_draw_vertices::Bool = false,
                                    f_matching::Bool = false
                                    )
    c,cr,bb = cairo_setup( filename, list, f_pdf );

    u_width::Float64 = 1024.0 * (BBox_width( bb) / 800.0);

    BBox_print( bb );
    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
#    set_line_width(cr, 10.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, u_width );
        if  len == 2  &&  count == 2
            set_source_rgb(cr, 0.0, 0.0, 1.0 );
        else
            set_source_rgb( cr, 0.0, 1.0, 0.0 );
        end

        draw_polygon( cr, poly );
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 1.5*u_width);
        set_source_rgb( cr, 1.0, 0.0, 1.0 );
        for  i  in 1:cardin( P )
            p = P[ i ];
            q = Q[ i ];

            move_to( cr, p[1], p[2] )
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 2.0*u_width);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, BBox_width( bb) / 200  );
        end
    end


    if  ( ! f_pdf )
        Cairo.write_to_png( c, filename );
    end
    Cairo.finish(c);
end


function  output_polygons_to_file_with_offsets(
    list::VecPolygon2F,
    loffs::VecVecFloat,
    filename,
    f_pdf::Bool,
    f_draw_vertices::Bool = false,
    f_matching::Bool = false
    )

    c,cr,bb = cairo_setup( filename, list, f_pdf );

    u_width::Float64 = 1024.0 * (BBox_width( bb) / 4500.0);

    BBox_print( bb );
    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
#    set_line_width(cr, 10.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    off::VecFloat = VecFloat();
    for  i in  eachindex(list)
        poly = list[ i ];
        f_off::Bool = false;
        if  ( i <= length( loffs ) )
            f_off = true;
            off = loffs[ i ];
        end
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, u_width );
        if  len == 2  &&  count == 2
            set_source_rgb(cr, 0.0, 0.0, 1.0 );
        else
            set_source_rgb( cr, 0.0, 1.0, 0.0 );
        end

        if  ( f_off )
            draw_polygon_w_offs( cr, poly, off );
        else
            draw_polygon( cr, poly );
        end
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 1.5*u_width);
        set_source_rgb( cr, 1.0, 0.0, 1.0 );
        for  i  in 1:cardin( P )
            p = P[ i ];
            q = Q[ i ];

            move_to( cr, p[1], p[2] )
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 2.0*u_width);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, BBox_width( bb) / 200  );
        end
    end


    if  ( ! f_pdf )
        Cairo.write_to_png( c, filename );
    end
    Cairo.finish(c);
end





#----------------------------------------------------------------
# Output the morphing to a pdf file
function  output_morphing( m::Morphing{N,T}, filename )  where {N,T}
    c,cr,bb = cairo_setup( filename, [ m.P, m.Q ], true );

    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_line_width(cr, 1.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    P = m.P;
    Q = m.Q;
    sol = m.sol;

    set_line_width(cr, 0.5);
    set_source_rgb(cr, 1.0, 0.0, 0.0 );

    nv::Int64 = cardin( sol );
    for  i in 1:nv
        loc::Point2I  = sol.pnts[ i ];

        po::Point{N,T} = P.pnts[ loc.x[ 1 ] ];
        qo::Point{N,T} = Q.pnts[ loc.x[ 2 ] ];

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


function  draw_frames( cr, sp::Segment2F, sq::Segment2F,
                       frames::Int64, P::Polygon2F, Q::Polygon2F, bb::BBox2F )

    delta::Float64 = 1.0 / (frames -1 );
    t::Float64 = 0.0;

    for i in 1:frames
        set_line_width(cr, 3.5);
        set_source_rgb(cr, 0.0, 0.8, 0.0 );
        draw_polygon( cr, P );
        set_source_rgb(cr, 0.0, 0.0, 1.0 );
        draw_polygon( cr, Q );


        set_line_width(cr, 10.5);
        set_source_rgb( cr, 1,0,0 );
        p::Point2F = Segment_get_on( sp, t );
        q::Point2F = Segment_get_on( sq, t );

        move_to( cr,  p[1], p[2]  )
        line_to( cr, q[1], q[2] );
        Cairo.stroke( cr );

        Cairo.show_page( cr );

        t = t + delta;
    end
end


function  compute_frames( pout, qout, total_frames )
#    println( "TOTAL FRAMES :", total_frames );
    lpout::Vector{Float64} = Polygon_prefix_lengths( pout )
    lqout::Vector{Float64} = Polygon_prefix_lengths( qout )

    lens::Vector{Float64} = lpout + lqout;

    steps::Int64 = length( lens ) - 1;
    total_length = last( lens );


    # compute how many frames for each leg...
    frames::Vector{Int64} = zeros( Int64, steps );
    # compute how many frames for each leg...
#    frames::Vector{Int64} = zeros( Int64, steps );
    acc = 0;
    for  i  in  1:steps
        leg_len = lens[ i + 1 ] - lens[ i ];
        if  ( leg_len <= 0 )
            continue;
        end
        acc += leg_len;
        num_pnts = round( Int64, total_frames * leg_len / total_length );
#        println( "num_pnts: ", num_pnts );
        if  num_pnts > 0
            acc = 0;
            frames[ i ] = max( num_pnts, 2 );
            continue;
        end

        # Very short edge...
        num_pnts = round( Int64, total_frames * acc / total_length );
        if  ( num_pnts == 0 )
            continue; # skip it....
        end
        acc = 0;
        frames[ i ] = num_pnts;
    end

    println( "Total frames #: ", sum( frames ) );
    return  frames, steps
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




mutable struct ContextMovie
    bb::BBox2F;
    bbo::BBox2F;
    iheight::Int64
    iwidth::Int64;
    frame_count::Int64;
    dir::String;
    f_show_vertices::Bool
end


mutable struct RecFrame
    frame::Int64
    p::Point2F
    q::Point2F
end
function   draw_image_record( cm::ContextMovie, P, Q, p, q,
                              vec::Vector{RecFrame} )
    cm.frame_count += 1;
    push!( vec, RecFrame( cm.frame_count, deepcopy( p ), deepcopy( q  ) ) );
end


"""
    draw_image_frame
"""
function   draw_image_frame( cm::ContextMovie, P, Q, rf::RecFrame )

    filename = cm.dir * "/" * @sprintf( "%06d.png", rf.frame )

    c = CairoRGBSurface(cm.iwidth, cm.iheight );
    cr = CairoContext(c);


    set_source_rgb(cr, 1, 1, 1);
    paint(cr);

    set_transform( cr, cm.iwidth, cm.iheight, cm.bbo );

    set_line_width(cr, 3.5);
    set_source_rgb(cr, 0.0, 0.8, 0.0 );
    draw_polygon( cr, P );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 1.0, 0.0 );
        draw_polygon_vertices( cr, P, BBox_width( cm.bb ) / 200 );
    end
    set_source_rgb(cr, 0.0, 0.0, 0.8 );
    draw_polygon( cr, Q );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 0.0, 1.0 );
        draw_polygon_vertices( cr, Q, BBox_width( cm.bb ) / 200 );
    end

    set_line_width(cr, 10.5);
    set_source_rgb( cr, 1,0,0 );

    move_to( cr,  rf.p[1], rf.p[2]  )
    line_to( cr, rf.q[1], rf.q[2] );
    Cairo.stroke( cr );

    print( "filename :", filename, "\r" );
    Cairo.write_to_png( c, filename );
    #Cairo.show_page( cr );
end



function  draw_frames_images( cm::ContextMovie, sp::Segment2F, sq::Segment2F,
    frames::Int64, P::Polygon2F, Q::Polygon2F, vec::Vector{RecFrame}  )

    delta::Float64 = 1.0 / (frames -1 );
    t::Float64 = 0.0;

    for i in 1:frames
        p::Point2F = Segment_get_on( sp, t );
        q::Point2F = Segment_get_on( sq, t );

        draw_image_record( cm, P, Q, p, q, vec )
#        cm::ContextMovie, P, Q, p, q,
#                               )

        #draw_image_frame( cm, P, Q, p, q );

        t = t + delta;
    end
end

function  frames_generate( cm::ContextMovie, P, Q, vec_rf )
    for rf in vec_rf
        draw_image_frame( cm::ContextMovie, P, Q, rf )
    end
    return  1
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


    println( "Splitting to threads... Sit back ;)" );

    println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)
    println( "" );
    #exit( -1 );
    # Then we generate the frames (but this can be parallelized!
#    for  i in eachindex(vec_rf)
#        draw_image_frame( cm::ContextMovie, m.P, m.Q, vec_rf[ i ] )
#    end


    ####
    # We now call ffmpeg to create the movie...
    # ffmpeg -r 10 -i temp/ do.mp4

#    cmd = ( " -r 10 -i " *   "
#          * filename );
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
#    println( "ffmpeg $options" );
    println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               cm.dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);

    println( "Rencoding with handbrake..." );

    # HandBrakeCLI -Z  -i movie.mp4  -o movie_2.mp4
    options_2 = [ "-Z", "Android 1080p30", "-i", tmp_filename,
                 "-o", filename ];
    output_2 = read(pipeline( `HandBrakeCLI $options_2`,
                            stderr="/tmp/errs_2.txt" ), String);
    rmx( tmp_filename );
    if  isfile( filename )
        println( "Created... ", filename );
    end
end


mutable struct ORecFrame
    frame::Int64
    t::Float64
    R::Polygon2F
end


function  frames_generate_o( cm::ContextMovie, P, Q, vec_rf )
    for rf in vec_rf
        draw_image_frame_o( cm::ContextMovie, P, Q, rf )
    end
    return  1
end


function   encode_to_mp4( dir, filename )
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
    println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);

    println( "\n\n\n\n\n\n" );
    println( "Rencoding with handbrake..." );
    println( filename );
    println( "\n\n\n\n\n\n" );

    # HandBrakeCLI -Z  -i movie.mp4  -o movie_2.mp4
    options_2 = [ "-Z", "Android 1080p30", "-i", tmp_filename,
                 "-o", filename ];
    output_2 = read(pipeline( `HandBrakeCLI $options_2`,
                            stderr="/tmp/errs_2.txt" ), String);
    rmx( tmp_filename );
    if  isfile( filename )
        println( "Created... ", filename );
    end
end


"""
    draw_image_frame_o
"""
function   draw_image_frame_o( cm::ContextMovie, P, Q, rf::ORecFrame )

    filename = cm.dir * "/" * @sprintf( "%06d.png", rf.frame )

    c = CairoRGBSurface(cm.iwidth, cm.iheight );
    cr = CairoContext(c);


    set_source_rgb(cr, 1, 1, 1);
    paint(cr);

    set_transform( cr, cm.iwidth, cm.iheight, cm.bbo );

    set_line_width(cr, 1);
    set_source_rgb(cr, 0.3, 0.3, 0.3 );
    draw_polygon( cr, P );
    draw_polygon( cr, Q );


    set_line_width(cr, 4);
    set_source_rgb(cr, 0.0, 0.8, 0.0 );
    draw_polygon( cr, rf.R );
    if   ( cm.f_show_vertices )
        set_source_rgb(cr, 0.0, 1.0, 0.0 );
        draw_polygon_vertices( cr, rf.R, BBox_width( cm.bb ) / 200 );
    end
    Cairo.stroke( cr );

    Cairo.write_to_png( c, filename );
end


function  frames_generate_o( cm::ContextMovie, P, Q,
                             vec_rf::Vector{ORecFrame} )
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

    println( "Splitting to threads... Sit back ;)" );

    println( Threads.nthreads() );

    chunks = Iterators.partition(vec_rf, length(vec_rf) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn frames_generate_o( cm, m.P, m.Q, chunk );
    end
    chunk_sums = fetch.(tasks)

    encode_to_mp4( cm.dir, filename );
end






# ????
# End of the frechet distance computation part
#################################################################
#################################################################
#################################################################



###################################################################3
# Compute an upper bound on the frechet distance between poly_a and
# its spine.
###################################################################3


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

#P::Polygon{N,T}, Q::Polygon{N,T}, Pe, Qe,
function  output_frechet_diagram( m::Morphing{N,T}, filename )  where {N,T}

    P_coords::Vector{T} = get_diagram_locs( m.pes, m.P );
    Q_coords::Vector{T} = get_diagram_locs( m.qes, m.Q );

    len = length( P_coords );
    poly = Polygon2F();
    for  i in 1:len
        Polygon_push_smart( poly, point( P_coords[ i ], Q_coords[ i ] ) );
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




function search_dir(path,key)
    #println( readdir( path ) );
    return  filter(x->occursin(key,x), readdir(path))
end




function  frechet_ve_r_mono_approx_dist( poly_a, poly_b, approx::Float64 )
    m = frechet_ve_r_mono_approx( poly_a, poly_b, approx )
    return  m.leash;
end

function  count_below_zero( arr )
    count = 0;
    for  x in arr
        if  x <= 0.0
            count = count + 1;
        end
    end
    return  count;
end

function   do_approx( poly_a, poly_b, approx::Float64 )
    @printf( "Approximating %3g", approx );
    m, ratio = frechet_c_approx( poly_a, poly_b, approx );
    println( "       : ", m.leash );
    return  m, ratio;
end

function  run_on_curves( poly_a, poly_b )
    println( " # vertices poly_a : ", cardin( poly_a ) );
    println( " # vertices poly_b : ", cardin( poly_b ) );

    do_approx( poly_a, poly_b, 10.0 );
    do_approx( poly_a, poly_b, 2.0 );
    do_approx( poly_a, poly_b, 1.1 );

    println( "Computing exact..." );
    m = frechet_c_compute( poly_a, poly_b );
    println( "Frechet distance        : ", m.leash );

    println( "Frechet discrete retractable..." );
    m = frechet_d_r_compute( poly_a, poly_b );

    println( "Frechet discrete... (quadratic time" );
    m = frechet_d_compute( poly_a, poly_b );
    println( "Computed!" );

    println( "Computing translated distance..." );
    P = deepcopy( poly_a );
    Polygon_translate!( P, -1.0 * first(poly_b) )
    m, iters = frechet_ve_r_compute( P, poly_b );
    println( "   Translated : ", m.leash, "  ", iters );
#    m, iters = frechet_ve_r_compute( poly_a, poly_b );
#    println( "Un Translated : ", m.leash, "  ", iters );

end

function  run_on_curves_files( filename_a, filename_b )
    println( "Reading curve: ", filename_a );
    poly_a = Polygon_read_plt_file( filename_a );
    println( "Reading curve: ", filename_b );
    poly_b = Polygon_read_plt_file( filename_b );

    return  run_on_curves( poly_a, poly_b );
end


function  do_example_1()
    poly_a,poly_b = example_1();
    output_polygons_to_flle( [poly_a, poly_b], "pdfs/curr_example.pdf", true );
    r_frechet_with_refinement( poly_a, poly_b, 40, true, false );
end

function  collect_arrows( xs, ys, vx, vy,
                          p::Point2F, q::Point2F )
    if  ( Dist(p, q ) == 0 )
        return;
    end
    push!( xs, p[1] );
    push!( ys, p[2] );

    push!( vx, q[1] - p[1] );
    push!( vy, q[2] - p[2] );
end

function  draw_arrow( plt, p::Point2F, q::Point2F, width = 1, ucolor=:pink )
    plot!(plt, [p[1],q[1]],[p[2],q[2]],
        color=ucolor,
        linewidth=width,
        label=:none, ticks=false, showaxis=false, grid=:false,
        legend=false, framestyle=:none,
        arrow=true
    )
end

"""
    diagram_get_ev_loc

Return the location in the Frechet free space diagram, of the edge-vertex
event between the edge P{i:i+1] and the vertex Q[j]
"""
function  diagram_get_ev_loc( P::Polygon{N,T}, Q::Polygon{N,T},
    pl::Vector{T}, ql::Vector{T},
    i::Int64, j::Int64 ) where  {N,T}
    seg = Segment( P[ i ], P[ i + 1 ] );
    p = Segment_nn_point( seg, Q[ j ] );
    x = pl[ i ] + Dist( P[ i ], p )
    y = ql[ j ]

    return point( x, y );
end

function  diagram_get_ve_loc( P::Polygon{N,T}, Q::Polygon{N,T},
    pl::Vector{T}, ql::Vector{T},
    i::Int64, j::Int64 ) where {N,T}
    seg = Segment( Q[ j ], Q[ j + 1 ] );
    q = Segment_nn_point( seg, P[ i ] );
    x = pl[ i ]
    y = ql[ j ] + Dist( Q[ j ], q )

    return point( x, y );
end


m = Vector{FrechetDist.Morphing2F}();







function  plot_curves_diagram( P::Polygon2F, Q::Polygon2F,
                               filename_diagram,
                               f_draw_c::Bool = false,
                               f_draw_ve::Bool = true,
                               f_draw_graph::Bool = true,
                               f_m_out_defined = false,
                               m_out = nothing,
                               title::String = "",
                               f_draw_graph_only::Bool = false
                               )

    println( "Getting ready to draw heatmap/graph/curves..." );

    len_P = Polygon_length( P )
    len_Q = Polygon_length( Q )

    pl = Polygon_prefix_lengths( P );
    ql = Polygon_prefix_lengths( Q );

    x_range = range(0, len_P, length = 200 )
    y_range = range(0, len_Q; length = 200 )

    function  fz(x,y)
        p = Polygon_get_point_on( P, pl, x );
        q = Polygon_get_point_on( Q, ql, y );

        return Dist( p, q );
    end

    println( "Computing heat map..." );
    ### f_draw_graph_only
    if  f_draw_graph_only
        plt = plot( x_range, y_range, 0,
                    left_margin = 0 * Plots.mm,
                    bottom_margin=0*Plots.mm,
                    right_margin = 0.02 * Plots.mm,
                    ticks = false, showaxis = false, framestyle=:none,
                    dpi = 200 );
    else
        plt = heatmap( x_range, y_range, fz,
                       color = :haline,
                       left_margin = 0 * Plots.mm,
                       bottom_margin=0*Plots.mm,
                       right_margin = 0.02 * Plots.mm,
                       ticks = false, showaxis = false, framestyle=:none,
                       dpi = 200 );
    end
    if  ( length( title ) > 0 )
        title!( plt, title );
    end
    println( "Heat map drawing done..." );

    function  draw_solution( plt )
        if  ( f_draw_c )
            m_c = frechet_c_compute( P, Q );
            p_c_diag = Morphing_extract_prm( m_c );
            m_c_diag = Polygon_as_matrix( p_c_diag );
            plot!(plt, m_c_diag[1,:], m_c_diag[2,:],
                linewidth=2, label=:none, ticks=false,
                showaxis=false, grid=:false,
                legend=false, framestyle=:none, lc=:red);
        end
        if  ( f_draw_ve )
            m_ve = frechet_ve_r_compute( P, Q );
            p_ve_diag = Morphing_extract_prm( m_ve );
            m_ve_diag = Polygon_as_matrix( p_ve_diag );
            plot!(plt, m_ve_diag[1,:], m_ve_diag[2,:],
                linewidth=4,
                label=:none, ticks=false, showaxis=false,
                grid=:false, legend=false, framestyle=:none, lc=:red);
        end
        println( "Drawing arrows..." );
    end

    function  draw_grid( plot )
        qlx = (ql[2:end-1])'
        plx = (pl[2:end-1])'

        if  f_draw_graph_only
            plot!( plt, [0; len_P], [qlx;qlx], lw=0.5, lc=:lightblue,
                   legend=false)
            plot!( plt,  [plx;plx], [0; len_Q], lw=0.5, lc=:lightblue,
                   legend=false)
        else
            plot!( plt, [0; len_P], [qlx;qlx], lw=2, lc=:black, legend=false)
            plot!( plt,  [plx;plx], [0; len_Q], lw=2, lc=:black, legend=false)
        end
    end
    function  draw_graph( plt )
        draw_solution( plt )

        p_last = cardin(P)-1;
        q_last = cardin(Q)-1;

        pnt_start::Point2F = point( 0.0, 0.0 );
        pnt_end::Point2F = point( Polygon_length( P ), Polygon_length( Q ) );

        pnts = Vector{Point2F}();
        s_a = point( 0.0, 0.0 );
        s_b = point( 0.0, 0.0 );
        t_a = point( 0.0, 0.0 );
        t_b = point( 0.0, 0.0 );
        counter::Int64 = 0;
        width = 1;
        ucolor = :pink;
        if   f_draw_graph_only
            width = 2;
            ucolor = :black;
        end

        xs = Float64[];
        ys = Float64[];
        vx = Float64[];
        vy = Float64[];
        for i in 1:p_last
            for j in 1:q_last
                s_a = diagram_get_ev_loc( P, Q, pl, ql, i, j )
                s_b = diagram_get_ve_loc( P, Q, pl, ql, i, j )
#                end

                f_use_s_b::Bool = ( i > 1 )  ||  ( j == 1 );

                f_use_t_a::Bool = (j < q_last);
                f_use_t_b::Bool = (i < p_last);
                if  ( i == p_last )  &&  ( j == q_last )
                    t_a = t_b = point( last( pl ), last( ql ) );
                    f_use_t_a = f_use_t_b = true;
                else
                    t_a = diagram_get_ev_loc( P, Q, pl, ql, i, j + 1 )
                    t_b = diagram_get_ve_loc( P, Q, pl, ql, i + 1, j  )
                end

                f_use_t_a && collect_arrows( xs, ys, vx, vy, s_a, t_a );
                f_use_t_b && collect_arrows( xs, ys, vx, vy, s_a, t_b );
                f_use_t_a && f_use_s_b && collect_arrows( xs, ys, vx, vy, s_b, t_a );
                f_use_t_b && f_use_s_b && collect_arrows( xs, ys, vx, vy, s_b, t_b );

                if  ( i == 1 ) &&  ( j == 1 )
                    collect_arrows( xs, ys, vx, vy, pnt_start, s_a );
                    collect_arrows( xs, ys, vx, vy, pnt_start, s_b );
                end
                #draw_arrow( plt, s_a, t_a, width, ucolor );
                #draw_arrow( plt, s_a, t_b, width, ucolor );
                #draw_arrow( plt, s_b, t_a, width, ucolor );
                #draw_arrow( plt, s_b, t_b, width, ucolor );

                push!( pnts, s_a )
                f_use_s_b  &&  push!( pnts, s_b )
                f_use_t_a  &&  push!( pnts, t_a )
                f_use_t_b  &&  push!( pnts, t_b );
            end
        end

        println( xs );
        quiver!(plt, xs, ys, quiver=(vx, vy),
                color=ucolor,
                linewidth=width,
                label=:none, ticks=false, showaxis=false, grid=:false,
                legend=false, framestyle=:none );

        sort!( pnts );
        unique!( pnts );

        m_pnts = VecPnts_as_matrix( pnts );
#=        println( "-------------------------------" );
        println( m_pnts[1,:] );
        println( "-------------------------------" );
        println( m_pnts[2,:] );
        println( "-------------------------------" );
        =#
        scatter!( plt, m_pnts[1,:], m_pnts[2,:], mc=:yellow, lc=:darkgreen,
            ms=4, ma=1.0 );

        end_pnts = Vector{Point2F}();
        push!( end_pnts, pnt_start, pnt_end );

        m_end_pnts = VecPnts_as_matrix( end_pnts );


        scatter!( plt, m_end_pnts[1,:], m_end_pnts[2,:], mc=:red,
                  ms=4, ma=3.5 );
        #    display( plt )
    end
    cardi::Int64 = cardin(P) + cardin( Q );
    if  ( cardi < 2000 )
        if f_draw_c  ||  f_draw_ve  ||  f_draw_graph
            draw_grid( plt );
        end
    end
    if  ( cardi < 2000 )  &&  f_draw_graph
        println( "Drawing the graph..." );
        draw_graph( plt );
    else
        draw_solution( plt )
    end

    if  ( f_m_out_defined )
        p_c_diag = Morphing_extract_prm( m_out );
        m_c_diag = Polygon_as_matrix( p_c_diag );
        plot!(plt, m_c_diag[1,:], m_c_diag[2,:],
              linewidth=2, label=:none, ticks=false,
              showaxis=false, grid=:false,
              legend=false, framestyle=:none, lc=:red);
    end

    println( "Saving heatmap/graph... ", filename_diagram );
    savefig( plt, filename_diagram );

    println( "Outputing the curves..." );
#    output_polygons_to_file(  [P, Q], filename_curves, true );
    println( "Created: " );
#    println( "   ", filename_curves );
    println( "   ", filename_diagram );
end

function  do_example( polys )
    poly_a = polys[ 1 ];
    poly_b = polys[ 2 ];

    pdf_output_polygons( [poly_a, poly_b], "pdfs/curr_example.pdf" );
    r_frechet_with_refinement( poly_a, poly_b, 40, false, false );
end

########################################################################
# Main....
########################################################################


function  create_movie( P::Polygon{N,T}, Q::Polygon{N,T},
    total_frames, filename, m ) where {N,T}
#    println( "About to compute the Frechet distance..." );
#    println( "Done!" );
    output_frechet_movie_mp4( m, filename, total_frames );
#    println( "Created  : ", filename );
    return  m;
end

function chop_it( s::String, c )
    while ( last( s) == c )
        s= chop(s,tail=1);
    end
    return  String( s );
end

function  create_demo( title::String,
                       prefix::String,
                       poly_a, poly_b,
                       f_draw_c::Bool = false,
                       f_draw_ve::Bool = true,
                       note::String = "",
                       f_refinements::Bool = false )
    if  ! isdir( prefix )
        mkdir( prefix );
    end
    cardi = cardin( poly_a ) + cardin( poly_b );
    total_frames = min( 50 * (cardin( poly_a ) + cardin( poly_b )), 800 );

    filename_curves = prefix*"curves.pdf";
    println( "Outputing the curves..." );
    output_polygons_to_file(  [poly_a, poly_b], filename_curves, true );

    options_svg = [ prefix*"curves.pdf", prefix*"curves.svg" ];
    output = read(pipeline( `pdf2svg $options_svg`, stderr="/tmp/errs.txt" ),
        String);

    println( "Created: " );
    println( "   ", filename_curves );

    local P, Q, m_d, m_d_r, m_ve_r, m_refinments;
    local m_d_dtw, m_adtw;
    f_computed_d::Bool = false;
    f_sampled_10::Bool = false;

    #println( "1.BOGI!\n\n\n\n" );

    #####################################################################
    # Computes distances
    #
    # m_c: Continuous monotone Frechet
    # m_d: Discrete Frecheet (potentially sampled)
    # m_d_r: Discrete restructured Frechet
    #####################################################################

#    fcei = FrechetCExtraInfo( Polygon2F(), Polygon2F(), Vector{Float64}(),
#                              Vector{Float64}(), false );

    m_c = frechet_c_compute( poly_a, poly_b, true )
    if  f_draw_ve
        m_ve_r = frechet_ve_r_compute( poly_a, poly_b );
    end


    f_adtw::Bool = false;
    if  ( cardi < 5000 )
        f_adtw = true;
        if  ( cardi < 100 )
            f_sampled_10 = true;
            P = Polygon_sample_uniformly( poly_a, 10*cardin( poly_a ) );
            Q = Polygon_sample_uniformly( poly_b, 10*cardin( poly_a ) );
        else
            P = poly_a;
            Q = poly_b;
        end
        m_d = frechet_d_compute( P, Q );
        m_d_dtw = DTW_d_compute( P, Q );
        m_d_r = frechet_d_r_compute( P, Q );
        f_computed_d = true;
    end

    if  f_adtw
        m_adtw = ADTW_compute( poly_a, poly_b );
    end

    local m_refinements::Vector{Morphing2F} = Vector{Morphing2F}();
    #println( "BOGI!\n\n\n\n" );
    if   f_refinements
        m_refinemenets = Vector{Morphing2F}();
        frechet_mono_via_refinement_ext( poly_a, poly_b, m_refinements, true,
                                         1.000000001
                                      );
        println( "m_refinements.len: ", length(m_refinements ) );
    end


    #####################################################################
    # Creating movies/diagrams/etc
    #####################################################################

    output_polygons_to_file(  [poly_a, poly_b], prefix * "curves.png", false );
    create_movie( poly_a, poly_b, total_frames, prefix*"f_c_movie.mp4", m_c );

    is_mkdir( prefix*"ortho/" );

    output_ortho_frechet_movie_mp4(  m_c, prefix*"ortho/c.mp4" );

    if  f_draw_ve
        create_movie( poly_a, poly_b, total_frames,
                      prefix*"f_ve_r_movie.mp4", m_ve_r );

        PU, QU = Morphing_as_polygons( m_ve_r );
        output_polygons_to_file(  [PU, QU], prefix*"ve_matching.pdf",
                                  true, true, true );
        #XXX
    end

    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.pdf",
        false, false, false
    );
    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.png",
        false, false, false
    );

    if   f_refinements
        dir =  prefix * "steps/";
        is_mkdir( dir );
        for  i in eachindex( m_refinements )
            mx = m_refinements[ i ];
            mm = Morphing_monotonize( mx );

            err::Float64 = 100.0 * ((mm.leash - mx.leash) / mx.leash)
            if  ( ( err > 1.0 )  ||  ( err == 0 ) )
                s= @sprintf( "%.2f", err)
            else
                digs::Int64 = 2 + convert(Int64, ceil(log10(1/err)) )
                s= @sprintf( "%-40.*f", digs, err)

                s = chop_it( s, ' ' );
#                s= chop_it( s, '0' );
            end

            title_frm = @sprintf( "Frame %02d   Monotonicity error: %s%%",
                              i, s )
            png_out = dir*@sprintf( "%06d.png", i );
            plot_curves_diagram( poly_a, poly_b,
                                 png_out,
                                 false, false, false, true,
                                 mx,
                                 title_frm
                                 );
        end

        println( "Generating gif..." );

        options = [ "-delay", "50", "-loop", "0", dir * "*.png",
                    prefix * "refinements.gif" ];
        outx = read(pipeline( `convert $options`, stderr="/tmp/errs.txt" ),
                    String );
    end

    f_graph_drawn::Bool = false;
    if  ( cardi < 100 )
        f_graph_drawn = true;
        plot_curves_diagram( poly_a, poly_b, prefix*"g_diagram.pdf",
            false, false, true
        );
        plot_curves_diagram( poly_a, poly_b, prefix*"g_diagram.png",
            false, false, true
                             );
        plot_curves_diagram( poly_a, poly_b, prefix*"graph_only.pdf",
                             false, false, true, false, nothing, "", true );
        plot_curves_diagram( poly_a, poly_b, prefix*"graph_only.png",
                             false, false, true, false, nothing, "", true );
    end
    plot_curves_diagram( poly_a, poly_b, prefix*"c_diagram.pdf",
        true, false, false
    );
    plot_curves_diagram( poly_a, poly_b, prefix*"c_diagram.png",
        true, false, false
                         );
    if  f_draw_ve
        plot_curves_diagram( poly_a, poly_b, prefix*"ve_r_diagram.pdf",
                             false, true, true
                             );
        plot_curves_diagram( poly_a, poly_b, prefix*"ve_r_diagram.png",
                             false, true, true
                             );
    end
    # println( "   ", prefix*"f_c_movie.mp4" );



    Polygon_write_to_file( poly_a, prefix * "poly_a.txt" );
    Polygon_write_to_file( poly_b, prefix * "poly_b.txt" );

    if  f_computed_d
        output_polygons_to_file(  [P, Q],
                                  prefix * "polygons_sampled.png", false,
                                  true  );
    end
    #########################################################
    # Html file...
    #########################################################

    fl = open( prefix * "index.html", "w" )

    println( "Writing file\n\n\n\n\n" );

    write( fl, "<head>\n"
           *"<meta charset=\"UTF-8\">"
           *"<TITLE>$prefix</TITLE>\n"
           *"<script type=\"text/x-mathjax-config\">\n"
           * "MathJax.Hub.Config({ tex2jax: "
           * "{inlineMath: [[\'\$\',\'\$\'], [\'\\(','\\)']]}"
           * "});\n"
           * "</script>\n"
           * "<script type=\"text/javascript\"\n"
           * "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js"
           * "config=TeX-AMS-MML_HTMLorMML\">\n"
           * "</script>\n"
           * "<meta charset=\"UTF-8\">\n"
           * "</head>" )
    write( fl, "<body>\n" );
    write( fl, "<h1>", title, "</h1>\n" );


    println( "Cardinality of both curves : ", cardi );



    # Image of curves

    write( fl, "<hr>\n\n" );
    write( fl, "<img src=\"curves.png\" />\n" )
    write( fl, "<hr>\n" )

    if  ( length( note ) > 0 )
        write( fl, "\n" )
        write( fl, "<!--- NOTE START --->\n" )
        write( fl, note )
        write( fl, "\n" )
        write( fl, "<!--- NOTE END --->\n" )
        write( fl, "<hr>\n" )
    end

    ############################################################
    # Table...
    row_a = [ "<a href=\"poly_a.txt\">P</a>"
              string( cardin( poly_a ) )
              string( Polygon_length( poly_a ) ) ]
    row_a_x = permutedims( row_a );

    row_b =  ["<a href=\"poly_b.txt\">Q</a>"
              string( cardin( poly_b ) )
              string( Polygon_length( poly_b ) )
              ];
    row_b_x = permutedims( row_b );

    data = vcat( row_a_x, row_b_x  );

    pretty_table(fl, data;
                 header = (["Curves", "# Vertices", "Length"]),
                 allow_html_in_cells = true,
                 backend = Val(:html) )
    write( fl, "\n<hr>\n" )

    if  f_draw_ve
        B=fill("", (2,2) )
        B[1,1] = "Fréchet";
        B[1,2] = string( m_c.leash );
        B[2,1] = "VE Fréchet";
        B[2,2] = string( m_ve_r.leash ) ;
        pretty_table(fl, B;
                     header = (["Distance",  "Value"]),
                     allow_html_in_cells = true,
                     backend = Val(:html) )
    else
        println( fl, "Fréchet distance: ", m_c.leash, "\n\n" );
    end
    write( fl, "\n<hr>\n" )



    ###########################################################
    # Movie
    write( fl, "\n\n<h2>Animation: "*FrechetStr*" morphing</h2>" );
    write( fl, "\n\n <video controls autoplay " );
    write( fl, "   src=\"f_c_movie.mp4\" type=\"video/mp4\" />\n" );
    write( fl, "</video>\n" );
    println( fl, "<p>\n"
             * "This is the animation of the morphing computed \n"
             * " that is both continuous and monotone.<p>\n" );

    write( fl, "\n\n\n<hr>\n" )


    ###########################################################
    # Movie: ve_retractable
    #        "f_ve_r_movie.mp4"

    if  f_draw_ve
        write( fl, "\n\n<h2>Animation: VE Retractable "
               * FrechetStr *"</h2>" );
        write( fl, "\n\n <video controls autoplay " );
        write( fl, "   src=\"f_ve_r_movie.mp4\" type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );
        println( fl, "<p>\n"
                 * "This is the animation of the VE retractable morphing."
                 * " It is potentially not monotone (but it is continuous."
                 * "<p>\n"
                 * "\n\n\n<hr>\n" );
    end

    #---------------------------------------------------------

    write( fl, "<h2>Free space diagram heatmap:</h2>" )
    write( fl, "<img src=\"diagram.png\">\n" );

    write( fl, "\n<hr>\n" )

    if  ( f_draw_ve  )
            write( fl, "<h2>VE-"*FrechetStr*" Retractable solution:</h2>" )
        write( fl, "<img src=\"ve_r_diagram.png\">\n" );
    end

    write( fl, "<h2>" * FrechetStr * " cont+monotone solution:</h2>" )

    #        println( fl, "What the animation shows.<br>\n" );
    write( fl, "<img src=\"c_diagram.png\">\n" );
    write( fl, "<hr>\n" );


    if  ( f_computed_d )
        output_frechet_movie_mp4( m_d_dtw, prefix*"d_dtw.mp4",
                                  400, true );
        output_frechet_movie_mp4( m_d, prefix*"discrete_frechet.mp4",
                                  400, true );
        output_frechet_movie_mp4( m_d_r, prefix*"discrete_r_frechet.mp4",
                                  400, true );
        println( fl, "<h1>Discrete " * FrechetStr * "</h1>\n" );
        if  ( f_sampled_10 )
            println( fl,
                     "Generated by sampling 10 points along each edge...<br>\n\n" );
        end
        println( fl, "<img src=\"polygons_sampled.png\">\n\n" );
        println( fl, "<p>\n" );
        println( fl, "The resulting morphing - extended to continuous:\n" );
        write( fl, "<p>\n\n <video controls autoplay " );
        write( fl, "   src=\"discrete_frechet.mp4\" "
               *"type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );

        println( fl, "<p>\n"
                 *"Specifically, to get a smooth animation, the \n"
                 *" leash " * "is shown as moving continuously, by \n"
                 * " interpolating between the discrete locations.\n"
                 * "<p>\n\n" );

        println( fl, "<hr>\n\n" );
        println( fl, "<h3>The discrete retractable version</h3>\n\n" );
        write( fl, "\n\n" );
        write( fl, "<video controls autoplay " );
        write( fl, "   src=\"discrete_r_frechet.mp4\" "
               * " type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );

        println( fl, "<hr>\n\n" );
        println( fl, "<h3>The discrete dynamic time warping</h3>\n\n" );
        write( fl, "\n\n" );
        write( fl, "<video controls autoplay " );
        write( fl, "   src=\"d_dtw.mp4\" "
               * " type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );



        println( fl, "<hr>" );
        println( fl, "P # vertices: ", cardin( P ), "<br>" );
        println( fl, "P # vertices: ", cardin( Q ), "<br>" );
        println( fl, "DFréchet iters         : ", m_d.iters, "<br>" );
        println( fl, "Retract DFréchet iters : ", m_d_r.iters, "<br>" );
    end

    if  ( f_adtw )
        output_frechet_movie_mp4( m_adtw, prefix*"adtw.mp4",
            400, true );
    end
    
    if  ( f_refinements )
        println( fl, "<hr>" * "\n" );
        println( fl, "<h2>Refinement for removing monotonicity</h2>\n" );

        println( fl, "By introducing vertices in the middle of "
                 * "parts of the curves that are being traversed in "
                 * "the wrong direction, one can refine the solution, "
                 * "till effectively reaching the optimal monotone "
                 * "solution. This process is demonstrated below. "
                 * "As one can see, the error is negligible after a few "
                 * "(say four) iterations. After that, it becomes a bit "
                 * "pointless. \n" );
        println( fl, "<br>We emphasize that a monotone morphing can \n"
                 * "always\n"
                 * " be extracted, by monotonizing the current \n"
                 * "solution. \n"
                 * "This is easy and fast to do, and is the error \n"
                 * "accounted for in the below graphics.<br>" );
        write( fl, "<img src=\"refinements.gif\"><br>\n\n\n" );
    end

    write( fl, "\n\n<h2>Animation: "*FrechetStr*" morphing as morphing</h2>" );
    write( fl, "\n\n <video controls autoplay " );
    write( fl, "   src=\"ortho/c.mp4\" type=\"video/mp4\" />\n" );
    write( fl, "</video>\n" );

    if  f_adtw
        println( fl, "<hr>" * "\n" );
        println( fl, "<h2>ADTW</h2>\n" );
        write( fl, "\n\n <video controls autoplay " );
        write( fl, "   src=\"adtw.mp4\" type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );
    end

    println( fl, "<hr>\n" );
    dt=now();
    date_str = Dates.format(dt, "yyyy-mm-dd HH:MM:SS")
    write( fl, date_str );


    write( fl, "</body>\n" );


    close( fl );

end


function  create_demo_files( title::String,
                             prefix::String,
                             f_a::String,
                             f_b::String,
                             f_draw_c::Bool = false,
                             f_draw_ve::Bool = true,
                             note::String = ""
                             )
    poly_a = Polygon_read_plt_file( f_a );
    poly_b = Polygon_read_plt_file( f_b );
    create_demo( title, prefix, poly_a, poly_b, f_draw_c, f_draw_ve, note );

end

gurl_base::String = "https://www.microsoft.com/en-us/research/publication/";
gurl_suffix::String = "geolife-gps-trajectory-dataset-user-guide/";

gurl::String = gurl_base * gurl_suffix;

gemb::String = "<a href=\"" * gurl * "\">GeoLife GPS Trajectories</a>";

function  gen_example_12()
    if  ( isfile( "data/041/trajectory/20090429225015.plt" ) )
        create_demo_files( "Example of close curves (GPS tracks)",
            "output/12/",
            "data/041/trajectory/20090429225015.plt",
            "data/041/trajectory/20090531225725.plt",
            true, false,
            "An example of two GPS tracks from " *
                gemb * " that are close together. \n" *
                 "This is an example where the retractable Fréchet\n" *
                 " algorithm axplores only tiny fraction of the diagam, \n" *
                 "yielding a near linear running time in this case.\n"
        );
    end
end


function  gen_example_6()
    poly_a,poly_b = example_6();
    create_demo( "Example 6: Refinement in action",
                 "output/06/",
                 poly_a,poly_b,
                 true, true,
                 "Zig-zag heavy example that shows the algorithm computing\n"
                 * " the exact continuous monotone " * FrechetStr
                 * " morphing, using refinement.\n",
                 true
                 );
end


function  gen_example_14()
    poly_a,poly_b = example_14();
    create_demo( "Example 14 the lone zigzag: Refinement in action",
                 "output/14/",
                 poly_a,poly_b,
                 true, true,
                 "Zig-zag heavy example that shows the algorithm computing\n"
                 * " the exact continuous monotone " * FrechetStr
                 * " morphing, using refinement.\n",
                 true
                 );
end

function  gen_example_15()
    poly_a,poly_b = example_15();
    create_demo( "Example 15",
                 "output/15/",
                 poly_a,poly_b,
                 true, true,
                 "Simple geometric shapes",
                 true
                 );
end


function  gen_example_16()
    poly_a,poly_b = example_16();
    create_demo( "Example 16", "output/16/", poly_a,poly_b,
                 true, true,
                 "Spirals",
                 true
                 );
end

function  gen_example_17()
    poly_a,poly_b = example_17_dtw();
    create_demo( "Example 17", "output/17/", poly_a,poly_b,
                 true, true,
                 "Example where DTW generates a different solution than"
                 *" discrete " * FrechetStr,
                 true
                 );
end

function  gen_example_1()
    poly_a,poly_b = example_1();
    create_demo( "Example 1", "output/01/",
                 poly_a,
                 poly_b,
                 false, true,
                 "A simple example demonstrating the main drawback \n"
                 * " of the regular discrete "*FrechetStr*" morphing \n"
                 * "which keeps the long leash after hitting the maximum \n"
                 * "length. The retractable discrete version on the other \n"
                 * "hand"
                 * " happily yields the \"correct\" result.\n"
                 );
end

function  is_rebuild( s::String )
    if  ( isdir( s ) )
        println( s, " already built." );
        return  false;
    end
    return  true
end

function  generate_examples()
    if  is_rebuild( "output/01" )
        gen_example_1()
    end
    if  is_rebuild( "output/02" )
        poly_a,poly_b = example_2();
        create_demo( "Example 2", "output/02/", poly_a,poly_b );
    end

    if  is_rebuild( "output/03" )
        poly_a,poly_b = example_3();
        create_demo( "Example 3", "output/03/", poly_a,poly_b );
    end

    if  is_rebuild( "output/04" )
        poly_a,poly_b = example_4();
        create_demo( "Example 4", "output/04/", poly_a,poly_b );
    end

    if  is_rebuild( "output/05" )
        poly_a,poly_b = example_5();
        create_demo( "Example 5", "output/05/", poly_a,poly_b );
    end

    if  is_rebuild( "output/06" )
        gen_example_6()
    end

    if  is_rebuild( "output/07" )
        poly_a,poly_b = example_7();
        create_demo( "Example 7", "output/07/", poly_a,poly_b );
    end

    if  is_rebuild( "output/08" )
        poly_a,poly_b = example_8_ext(3);
        create_demo( "Example 8 (3)", "output/08/", poly_a,poly_b );
    end

    if  is_rebuild( "output/09" )
        poly_a,poly_b = example_9(  3, 4 );
        create_demo( "Example 9 (double zig-zag 3+4)", "output/09/",
                     poly_a,poly_b );
    end

    if  is_rebuild( "output/10" )
        poly_a,poly_b = example_10( 3, 4);
        create_demo( "Example 10", "output/10/", poly_a,poly_b );
    end

    if   is_rebuild( "output/11" )
        if  ( isfile( "data/010/trajectory/20080928160000.plt" ) )
            create_demo_files( "Example of big curves (GPS tracks)",
                               "output/11/",
                               "data/010/trajectory/20080928160000.plt",
                               "data/010/trajectory/20081219114010.plt",
                               true, false );
        end
    end

    if  is_rebuild( "output/12" )
        gen_example_12();
    end

    if  is_rebuild( "output/13" )
        poly_a,poly_b = example_13();
        create_demo( "Example 13 (7 fixed)", "output/13/", poly_a,poly_b );
    end

    if  is_rebuild( "output/14" )
        gen_example_14()
    end

    if  is_rebuild( "output/15" )
        gen_example_15()
    end

    if  is_rebuild( "output/16" )
        gen_example_16()
    end
    if  is_rebuild( "output/17" )
        gen_example_17()
    end
end


if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );
if   num_args == 0
#    gen_example_6();

    generate_examples();
    exit( 0 );
end

if   num_args == 2
    create_demo_files( "output/test/", "output/test/",
                       ARGS[1], ARGS[2], true, false );
    exit( 0 );
end

if   num_args == 3
    create_demo_files( "empty_title", ARGS[1], ARGS[2], ARGS[3], true, false );
    exit( 0 );
end


exit( -1 );

###############################################################################
###############################################################################
###############################################################################

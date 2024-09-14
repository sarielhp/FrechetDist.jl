# Graphics.jl

using Parameters
using Distributions;
using Cairo
using LinearAlgebra
using Printf
using Plots
using PrettyTables
using Dates


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


function  draw_segments( cr, segs )
    for  s in segs
        p = s.p;
        q = s.q;
        move_to(cr, p[ 1 ], p[ 2 ] );
        line_to(cr, q[ 1 ], q[ 2 ] );
        Cairo.stroke( cr );
    end
end


function  draw_polygon( cr, P, f_close::Bool = false )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        p = P.pnts[ i - 1 ];
        #p = po.x;
        q = P.pnts[ i ];
        #q = qo.x;
        move_to( cr,  p[1], p[2] )
        line_to( cr, q[1], q[2] );
        #println( q[1], " ", q[2] );
    end
    if  (nv > 0 )  &&  ( f_close )
        q = P.pnts[ 1 ];
        line_to( cr, q[1], q[2] );
    end
    Cairo.stroke( cr );
end


function  draw_points( cr, P::Vector{Point2F}, rad::Float64 = 0.005 )
    nv::Int64 = length( P );
    #println( "nv:", nv );
    for  i in 1:nv
        p = P[ i ];
        move_to( cr, p[1] - rad , p[2] )
        
        arc( cr, p[1], p[2], rad, 0.0, 2 * 3.1415 );
        stroke_preserve(cr);
        fill(cr);
#        move_to( cr,  p[1], p[2] )
#        line_to( cr,  p[1]+0.01, p[2] )
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

    println( BBox_width( bb) );
    u_width::Float64 = 3.0; #1024.0 * (BBox_width( bb) / 200.0);
#    u_width = 3;
#    exit( -1 );

    #BBox_print( bb );
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

    #println( "Total frames #: ", sum( frames ) );
    return  frames, steps
end


<<<<<<< HEAD
=======
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
    frames, steps = compute_frames( pout, qout, total_frames )

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

>>>>>>> 7e19abc (minor changes)



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

    #print( "filename :", filename, "\r" );
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





function   encode_to_mp4( dir, filename )
    tmp_filename = "tmp.mp4";
    rmx( filename );
    rmx( tmp_filename );
    #println( "Encoding movie with ffmpeg..." );
    options = [ "-r", "10", "-i",
               dir * "/" * "%06d.png",
               "-c:v", "libx264", tmp_filename ];
    output = read(pipeline( `ffmpeg $options`, stderr="/tmp/errs.txt" ),
        String);


    #println( "\n\n\n\n\n\n" );
    #println( "Rencoding with handbrake..." );
    #println( filename );
    #println( "\n\n\n\n\n\n" );

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

mutable struct ORecFrame
    frame::Int64
    t::Float64
    R::Polygon2F
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








# ????
# End of the frechet distance computation part
#################################################################
#################################################################
#################################################################

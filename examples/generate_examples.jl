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
    f_draw_vertices::Bool = false )
    c,cr,bb = cairo_setup( filename, list, f_pdf );

    BBox_print( bb );
    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_line_width(cr, 10.0);
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        #println( count, " ", len );
        set_line_width(cr, 2.00);
        if  len == 2  &&  count == 2
            set_source_rgb(cr, 0.0, 0.0, 1.0 );
        else
            set_source_rgb( cr, 0.0, 1.0, 0.0 );
        end

        draw_polygon( cr, poly );
        if  ( f_draw_vertices )
            set_line_width(cr, 8.00);
            set_source_rgb( cr, 1.0, 0.0, 0.0 );
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

function  draw_arrow( plt, p::Point2F, q::Point2F )
    plot!(plt, [p[1],q[1]],[p[2],q[2]],
        color=:pink,
        linewidth=1,
        label=:none, ticks=false, showaxis=false, grid=:false,
        legend=false, framestyle=:none,
        arrow=false
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

function  plot_curves_diagram( P, Q,
    filename_diagram,
    f_draw_c::Bool = false,
    f_draw_ve::Bool = true,
    f_draw_graph::Bool = true )

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
    plt = heatmap( x_range, y_range, fz, color = :haline,#:copper, #:thermal,
                   left_margin = 0 * Plots.mm,
                   bottom_margin=0*Plots.mm,
                   right_margin = 0.02 * Plots.mm,
        ticks = false, showaxis = false, framestyle=:none)
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

        plot!( plt, [0; len_P], [qlx;qlx], lw=2, lc=:black, legend=false)
        plot!( plt,  [plx;plx], [0; len_Q], lw=2, lc=:black, legend=false)
    end
    function  draw_graph( plt )
        draw_solution( plt )

        p_last = cardin(P)-1;
        q_last = cardin(Q)-1;
        pnts = Vector{Point2F}();
        s_a = point( 0.0, 0.0 );
        s_b = point( 0.0, 0.0 );
        t_a = point( 0.0, 0.0 );
        t_b = point( 0.0, 0.0 );
        counter::Int64 = 0;
        for i in 1:p_last
            for j in 1:q_last
                if  ( i == 0 ) &&  ( j == 0 )
                    s_a = s_b = point( 0.0, 0.0 );
                else
                    s_a = diagram_get_ev_loc( P, Q, pl, ql, i, j )
                    s_b = diagram_get_ve_loc( P, Q, pl, ql, i, j )
                end
                if  ( i == p_last )  &&  ( j == q_last )
                    t_a = t_b = point( last( pl ), last( ql ) );
                else
                    t_a = diagram_get_ev_loc( P, Q, pl, ql, i, j + 1 )
                    t_b = diagram_get_ve_loc( P, Q, pl, ql, i + 1, j  )
                end
                draw_arrow( plt, s_a, t_a );
                draw_arrow( plt, s_a, t_b );
                draw_arrow( plt, s_b, t_a );
                draw_arrow( plt, s_b, t_b );

                push!( pnts, s_a, s_b, t_a, t_b );
            end
        end
        sort!( pnts );
        unique!( pnts );

        m_pnts = VecPnts_as_matrix( pnts );
#=        println( "-------------------------------" );
        println( m_pnts[1,:] );
        println( "-------------------------------" );
        println( m_pnts[2,:] );
        println( "-------------------------------" );
        =#
        scatter!( plt, m_pnts[1,:], m_pnts[2,:], mc=:green, ms=2, ma=0.5 );

        # Copy the first and last point of m to u...
        u = hcat( m_pnts[:,1], m_pnts[ :, size( m_pnts, 2 ) ] );

        scatter!( plt, u[1,:], u[2,:], mc=:red, ms=4, ma=3.5 );
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




function  create_demo( title::String, prefix, poly_a, poly_b,
                       f_draw_c::Bool = false,
                       f_draw_ve::Bool = true,
                       note::String = "" )
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

    local P, Q, m_d, m_d_r, m_ve_r
    f_computed_d::Bool = false;
    f_sampled_10::Bool = false;

    #####################################################################
    # Computes distances
    #
    # m_c: Continuous monotone Frechet
    # m_d: Discrete Frecheet (potentially sampled)
    # m_d_r: Discrete restructured Frechet
    #####################################################################
    m_c = frechet_c_compute( poly_a, poly_b )
    if  f_draw_ve
        m_ve_r = frechet_ve_r_compute( poly_a, poly_b );
    end

    if  ( cardi < 5000 )
        if  ( cardi < 100 )
            f_sampled_10 = true;
            P = Polygon_sample_uniformly( poly_a, 10*cardin( poly_a ) );
            Q = Polygon_sample_uniformly( poly_b, 10*cardin( poly_a ) );
        else
            P = poly_a;
            Q = poly_b;
        end
        m_d = frechet_d_compute( P, Q );
        m_d_r = frechet_d_r_compute( P, Q );
        f_computed_d = true;
    end



    #####################################################################
    # Creating movies/diagrams/etc
    #####################################################################

    output_polygons_to_file(  [poly_a, poly_b], prefix * "curves.png", false );
    create_movie( poly_a, poly_b, total_frames, prefix*"f_c_movie.mp4", m_c );

    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.pdf",
        false, false, false
    );
    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.png",
        false, false, false
    );

    f_graph_drawn::Bool = false;
    if  ( cardi < 100 )
        f_graph_drawn = true;
        plot_curves_diagram( poly_a, poly_b, prefix*"g_diagram.pdf",
            false, false, true
        );
        plot_curves_diagram( poly_a, poly_b, prefix*"g_diagram.png",
            false, false, true
        );
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
    println( "   ", prefix*"f_c_movie.mp4" );



    Polygon_write_to_file( poly_a, prefix * "poly_a.txt" );
    Polygon_write_to_file( poly_a, prefix * "poly_b.txt" );

    if  f_computed_d
        output_polygons_to_file(  [P, Q],
                                  prefix * "polygons_sampled.png", false,
                                  true  );
    end
    #########################################################
    # Html file...
    #########################################################

    open( prefix * "index.html", "w" ) do fl
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
            write( fl, note )
            write( fl, "\n" )
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
        write( fl, "\n\n<h2>Animation of the Frechet morphing</h2>" );
        write( fl, "\n\n <video controls autoplay " );
        write( fl, "   src=\"f_c_movie.mp4\" type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );
        println( fl, "<p>\n"
        * "This is the anomation of the morphing computed \n"
        * " that is both continuous and monotone.<p>\n" );


        write( fl, "\n\n\n<hr>\n" )


        write( fl, "<h2>Free space diagram heatmap:</h2>" )
        write( fl, "<img src=\"diagram.png\">\n" );

        write( fl, "\n<hr>\n" )

        if  ( f_draw_ve  )
            write( fl, "<h2>VE-Frechet Retractable solution:</h2>" )
            write( fl, "<img src=\"ve_r_diagram.png\">\n" );
        end

        write( fl, "<h2>Frechet cont+monotone solution:</h2>" )

#        println( fl, "What the animation shows.<br>\n" );
        write( fl, "<img src=\"c_diagram.png\">\n" );
        write( fl, "<hr>\n" );


        if  ( f_computed_d )
            output_frechet_movie_mp4( m_d, prefix*"discrete_frechet.mp4",
                400, true );
            output_frechet_movie_mp4( m_d_r, prefix*"discrete_r_frechet.mp4",
                400, true );
            println( fl, "<h1>Discrete Frechet</h1>\n" );
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
                *"Specifically, to get a smooth animatino, the \n"
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

            println( fl, "<hr>" );
            println( fl, "P # vertices: ", cardin( P ), "<br>" );
            println( fl, "P # vertices: ", cardin( Q ), "<br>" );
            println( fl, "DFréchet iters         : ", m_d.iters, "<br>" );
            println( fl, "Retract DFréchet iters : ", m_d_r.iters, "<br>" );
        end



        println( fl, "<hr>\n" );
        dt=now();
        date_str = Dates.format(dt, "yyyy-mm-dd HH:MM:SS")
        write( fl, date_str );


        write( fl, "</body>\n" );

    end


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
    create_demo( title, prefix, poly_a, poly_b, f_draw_c, f_draw_ve );

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



function  generate_examples()
    poly_a,poly_b = example_1();
    create_demo( "Example 1", "output/01/", poly_a,poly_b );

    poly_a,poly_b = example_2();
    create_demo( "Example 2", "output/02/", poly_a,poly_b );

    poly_a,poly_b = example_3();
    create_demo( "Example 3", "output/03/", poly_a,poly_b );

    poly_a,poly_b = example_4();
    create_demo( "Example 4", "output/04/", poly_a,poly_b );

    poly_a,poly_b = example_5();
    create_demo( "Example 5", "output/05/", poly_a,poly_b );

    poly_a,poly_b = example_6();
    create_demo( "Example 6", "output/06/", poly_a,poly_b );
    poly_a,poly_b = example_7();
    create_demo( "Example 7", "output/07/", poly_a,poly_b );
    poly_a,poly_b = example_8_ext(3);
    create_demo( "Example 8 (3)", "output/08/", poly_a,poly_b );

    poly_a,poly_b = example_9(  3, 4 );
    create_demo( "Example 9 (double zig-zag 3+4)", "output/09/",
                 poly_a,poly_b );

    poly_a,poly_b = example_10( 3, 4);
    create_demo( "Example 10", "output/10/", poly_a,poly_b );

    if  ( isfile( "data/010/trajectory/20080928160000.plt" ) )
        create_demo_files( "Example of big curves (GPS tracks)",
            "output/11/",
            "data/010/trajectory/20080928160000.plt",
            "data/010/trajectory/20081219114010.plt",
            true, false );
    end
    gen_example_12();
end

if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );
if   num_args == 0
    gen_example_12();

#    generate_examples();
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

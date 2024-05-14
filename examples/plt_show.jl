#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using Cairo
using FrechetDist
using FrechetDist.cg

VecFloat = Vector{Float64};
VecVecFloat = Vector{VecFloat};


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

    u_width::Float64 = 1024.0 * (BBox_width( bb) / 100.0);

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



function  plt_show( ARGS )
    num_args = length( ARGS );

    list = VecPolygon2F();
    for  i in 1:num_args
        poly_a = Polygon_read_plt_file( ARGS[ i ] );
        push!( list, poly_a );
    end
    output_polygons_to_file( list, "curves.pdf", true );
    println( "Generated curves.pdf" );
end

####################################################################

if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );


if   num_args == 0
    println( "plt_show [plt_file] ... [plt_file]" );
    exit( -1 );
end

plt_show( ARGS );

###############################################################################
###############################################################################
###############################################################################

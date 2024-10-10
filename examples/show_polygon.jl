#! /usr/bin/julia

push!(LOAD_PATH, pwd()*"/src/")

using Cairo
using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon
using Parameters

include( "graphics.jl" )

VecFloat = Vector{Float64};
VecVecFloat = Vector{VecFloat};


mutable struct  Affine2d
    x::Float64;
    y::Float64;
    scale::Float64;

    function  Affine2d()
        return  new( 0.0, 0.0, 1.0 );
    end
    function  Affine2d( _x, _y, _scale)
        return  new( _x, _y, _scale );
    end
end;

function  Base.:*(t::Affine2d, p::Point{D,T}) where  {D,T}
    return  npoint( ( p[ 1 ] + t.x)*t.scale,
        ( p[ 2 ] + t.y)*t.scale );
end


function  draw_polygon( cr, P, T::Affine2d )
    nv::Int64 = cardin( P );
    for  i in 2:nv
        #println( T * P[ i ] );
        po = T * P.pnts[ i - 1 ];
        qo = T * P.pnts[ i ];
        move_to( cr,  po[1], po[2] )
        line_to( cr, qo[1], qo[2] );
        #println( q[1], " ", q[2] );
    end
    Cairo.stroke( cr );
end

function  draw_polygon_vertices( cr, P, r::Float64 )
    nv::Int64 = cardin( P );
    #r = total_length( P ) / (100.0*cardin( P ));
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

    pb = npoint( pc[1], pa[2] );
    pd = npoint( pa[1], pc[2] )

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

    println( "--------------------" );
    BBox_print( bbo );
    println( "--------------------" );
    BBox_print( bb );
    println( "--------------------" );

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

    println( "Image dimensions: (", iwidth, ", ", iheight, ")" );
    return  iheight,iwidth;
end




function  set_transform( cr, iwidth::Int64, iheight::Int64,
                         bbo::BBox2F )
    xcal = convert( Float64, iwidth) / BBox_width( bbo, 1 );

    println( "Scaling: ", xcal );
    #Cairo.scale( cr, xcal, xcal );
    bl = BBox_bottom_left( bbo );
    println( "bl :", bl );
    #Cairo.translate( cr, -bl[ 1 ], -bl[ 2 ]);

    println( "IWIDTH: ", iwidth );
    println( "TRANSLATION: ", bl[1], " , ", bl[2 ] );
    return  Affine2d( -bl[ 1 ], -bl[ 2 ], xcal );
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

    T = set_transform( cr, iwidth, iheight, bbo );

    if  ( ! f_pdf )
        set_source_rgb(cr, 1, 1, 1);
        paint(cr);
    end

    return  c,cr,bb, T;
end



function  output_polygons_to_file(  list::VecPolygon2F, filename,
                                    f_pdf::Bool,
                                    f_draw_vertices::Bool = false,
                                    f_matching::Bool = false
                                    )
    c,cr,bb, T = cairo_setup( filename, list, f_pdf );


    u_width::Float64 = 1024.0 * (BBox_width( bb) / 100.0);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        println( count, "/", len, "   |P|:", cardin( poly ) );

        set_line_width(cr, 2.0 + 3.0 / count );
        #println( "IND = ", ind );
        clr = get_color_rgb( count );

        #=
        if  ( count == 3 ) 
            set_line_width(cr, 1.0 );
            clr = get_color_rgb( 2 );
            set_dash( cr,  [4.0,4.0, 4.0]);
        end
        =#
        set_source_rgb( cr, clr... );
        draw_polygon( cr, poly, T );
    end

    if  ( f_matching )  &&  ( cardin( list[ 1 ] ) ==  cardin( list[ 2 ] ) )
        #println( "BOGI\n\n\n\n" );
        P = list[ 1 ];
        Q = list[ 2 ];
        set_line_width(cr, 0.3*u_width);
        set_source_rgb( cr, 0.0, 0.5, 0.0 );
        for  i  in 1:cardin( P )
            p = T*P[ i ];
            q = T*Q[ i ];

            move_to( cr, p[1], p[2] )
            set_dash( cr,  [4.0, 21.0, 2.0]);
            line_to( cr, q[1], q[2] );
            Cairo.stroke( cr );
        end
    end

    if  ( f_draw_vertices )
        set_line_width(cr, 2.0);
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
    bb = BBox2F();
    for  i in 1:num_args
        println( "Reading: ", ARGS[ i ] );
        poly_a = read_file( ARGS[ i ] );
        push!( list, poly_a );
        BBox_bound( bb, poly_a );
    end
    p = BBox_bottom_left( bb );
    for  poly  in list
        Polygon_translate!( poly, p );
    end

    #for  poly  in list
    #    println( poly );
    #end

    output_polygons_to_file( list, "curves.pdf", true );
    println( "Generated curves.pdf" );
    output_polygons_to_file( list, "curves.png", false );
    println( "Generated curves.png" );

    #=
    output_polygons_to_file( list, "curves_m.pdf", true, true, true );
    println( "Generated curves.pdf" );
    output_polygons_to_file( list, "curves_m.png", false, true, true );
    println( "Generated curves.png" );
    =#
end

####################################################################

if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );


if   num_args == 0
    println( "show_polygon [txt/plt file] ... [txt/plt file]" );
    exit( -1 );
end

plt_show( ARGS );

############################################################################
############################################################################
############################################################################

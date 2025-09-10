#! julia

push!(LOAD_PATH, pwd()*"/src/")

using Cairo
using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point
using FrechetDist.cg.polygon
using FrechetDist.cg.bbox
using Parameters

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
    pa = bottom_left( bb );
    pc = top_right( bb );

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

    bound( bb, list );
    expand!( bb, 1.05 );
    bbo = deepcopy( bb );
    println( typeof( bbo ) );
    expand!( bbo, 1.05 );

    println( "--------------------" );
    bbox.print( bbo );
    println( "--------------------" );
    bbox.print( bb );
    println( "--------------------" );

    return  bb, bbo
end

function  get_image_dims( bbo )

    width::Float64 = 1024.0;
    theight::Float64 = 0.0;

    while ( true )
        theight = width * bbox.width( bbo, 2 ) / bbox.width( bbo, 1 );
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
    xcal = convert( Float64, iwidth) / bbox.width( bbo, 1 );

    println( "Scaling: ", xcal );
    #Cairo.scale( cr, xcal, xcal );
    bl = bottom_left( bbo );
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


function  output_polygons_hierarchy_to_file(  list::VecPolygon2F, filename,
    f_pdf::Bool,
    f_draw_vertices::Bool = false,
    f_matching::Bool = false
                                    )
    c,cr,bb, T = cairo_setup( filename, list, f_pdf );

    u_width::Float64 = 1024.0 * (bbox.width( bb) / 100.0);

    set_source_rgb(cr,0.9,0.9,0.9);    # light gray
    set_source_rgba(cr, 1, 0.2, 0.2, 0.6);

    len = length( list );
    count::Int64 = 0;
    for  poly in  list
        count = count + 1;
        println( count, "/", len, "   |P|:", cardin( poly ) );

        palette = frechet_palette( poly );

        limit = 5;
        len = cardin( poly );
        lg = round( Int64, log( len ) / log( 2.0 ) );
        println( "LG: ", lg );
        delta = round( Int64, lg / limit );
        lvl::Int64 = lg;
        println( "\n\n\n" );
        for  i in  limit:-1:0
            set_line_width(cr, 25.0/ ( 1 +  limit - i + 1) );
            ( i == 0 ) &&  set_source_rgb(cr, 0.0, 0.0, 1.0 );
            ( i == 1 ) &&  set_source_rgb(cr, 1.0, 0.0, 0.0 );
            ( i == 2 ) &&  set_source_rgb(cr, 0.0, 1.0, 0.0 );
            ( i == 3 ) &&  set_source_rgb(cr, 0.0, 1.0, 1.0 );
            ( i == 4 ) &&  set_source_rgb(cr, 1.0, 1.0, 0.0 );
            ( i == 4 ) &&  set_source_rgb(cr, 1.0, 0.0, 1.0 );

            println( "lvl: ", lvl );
            Q = frechet_palette_level( poly, palette, lvl )[1];
            println( cardin( Q ) );
            draw_polygon( cr, Q, T );
            lvl = max( lvl - delta, 0 );
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
        set_line_width(cr, 2.0);
        set_source_rgb( cr, 1.0, 0.0, 0.0 );
        for  poly in  list
            draw_polygon_vertices( cr, poly, bbox.width( bb) / 200  );
        end
    end


    if  ( ! f_pdf )
        Cairo.write_to_png( c, filename );
    end
    Cairo.finish(c);
end



function  plt_show_ph( ARGS )
    num_args = length( ARGS );

    list = VecPolygon2F();
    bb = BBox2F();
    for  i in 1:num_args
        println( "Reading: ", ARGS[ i ] );
        poly_a = polygon.read_file( ARGS[ i ] );
        push!( list, poly_a );
        bound( bb, poly_a );
    end
    p = bottom_left( bb );
    for  poly  in list
        Polygon_translate!( poly, p );
    end

    #for  poly  in list
    #    println( poly );
    #end

    output_polygons_hierarchy_to_file( list, "curves.pdf", true );
    println( "Generated curves.pdf" );
end

####################################################################

if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );


if   num_args == 0
    println( "show_polygon_hierarchy.jl [txt/plt file]" );
    exit( -1 );
end

plt_show_ph( ARGS );

############################################################################
############################################################################
############################################################################

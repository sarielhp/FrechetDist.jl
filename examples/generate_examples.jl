#! /bin/julia

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
using FrechetDist.cg.point
using FrechetDist.cg.segment
using FrechetDist.cg.polygon
#using cg

include( "fr_examples.jl" )

include( "graphics.jl" )

FrechetStr::String = "Fréchet";


function  test_command( cmds... )
    for  cmd  in cmds
#        println( "cmd: ", cmd );
        if  ( Sys.which( cmd ) == nothing )
            println( "Error:\n\t"*
                "Required system commad [", cmd,"] not found.\n\n"*
                 "Please install, and try again!\n\n"*
                "If you are missing \"convert\", this is provided "*
                    "by ImageMagic.\n\n" );
            exit( -1 );
        end
    end
end



function  is_mkdir( dir )
    if  ! isdir( dir )
        mkdir( dir );
    end
end


VecFloat = Vector{Float64};
VecVecFloat = Vector{VecFloat};





###################################################################3
# Compute an upper bound on the frechet distance between poly_a and
# its spine.
###################################################################3






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
    poly_a = polygon.read_file( filename_a );
    println( "Reading curve: ", filename_b );
    poly_b = polygon.read_file( filename_b );

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
    p = nn_point( seg, Q[ j ] );
    x = pl[ i ] + Dist( P[ i ], p )
    y = ql[ j ]

    return npoint( x, y );
end

function  diagram_get_ve_loc( P::Polygon{N,T}, Q::Polygon{N,T},
    pl::Vector{T}, ql::Vector{T},
    i::Int64, j::Int64 ) where {N,T}
    seg = Segment( Q[ j ], Q[ j + 1 ] );
    q = nn_point( seg, P[ i ] );
    x = pl[ i ]
    y = ql[ j ] + Dist( Q[ j ], q )

    return  npoint( x, y );
end


#m = Vector{FrechetDist.Morphing2F}();





function  plot_curves_diagram(
    P::Polygon2F, Q::Polygon2F,
    filename_diagram,
    f_draw_c::Bool = false,
    f_draw_ve::Bool = true,
    f_draw_graph::Bool = true,
    f_m_out_defined = false,
    m_out = nothing,
    title::String = "",
    f_draw_graph_only::Bool = false;
    f_draw_grid::Bool = false,
    f_draw_monotone = false,
    f_3d = false
)
###----------------------------------------------------------------------
### sub-function draw_graph start
    function  draw_graph( plt )
        draw_solution( plt )

        p_last = cardin(P)-1;
        q_last = cardin(Q)-1;

        pnt_start::Point2F = npoint( 0.0, 0.0 );
        pnt_end::Point2F = npoint( total_length( P ), total_length( Q ) );

        pnts = Vector{Point2F}();
        s_a = npoint( 0.0, 0.0 );
        s_b = npoint( 0.0, 0.0 );
        t_a = npoint( 0.0, 0.0 );
        t_b = npoint( 0.0, 0.0 );
        counter::Int64 = 0;
        width = 1;
        ucolor = :pink;
        if   f_draw_graph_only  &&  ( ! f_draw_grid )
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

                f_use_s_b::Bool = ( i > 1 )  ||  ( j == 1 );

                f_use_t_a::Bool = (j < q_last);
                f_use_t_b::Bool = (i < p_last);
                if  ( i == p_last )  &&  ( j == q_last )
                    t_a = t_b = npoint( last( pl ), last( ql ) );
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

                push!( pnts, s_a )
                f_use_s_b  &&  push!( pnts, s_b )
                f_use_t_a  &&  push!( pnts, t_a )
                f_use_t_b  &&  push!( pnts, t_b );
            end
        end

        f_debug && println( xs );
        quiver!(plt, xs, ys, quiver=(vx, vy),
            color=ucolor,
            linewidth=width,
            label=:none, ticks=false, showaxis=false, grid=:false,
            legend=false, framestyle=:none );

        sort!( pnts );
        unique!( pnts );

        m_pnts = VecPnts_as_matrix( pnts );
        scatter!( plt, m_pnts[1,:], m_pnts[2,:], mc=:yellow, lc=:darkgreen,
            ms=4, ma=1.0 );

        end_pnts = Vector{Point2F}();
        push!( end_pnts, pnt_start, pnt_end );

        m_end_pnts = VecPnts_as_matrix( end_pnts );


        scatter!( plt, m_end_pnts[1,:], m_end_pnts[2,:], mc=:red,
            ms=4, ma=3.5 );
        #    display( plt )
    end
### sub-function draw_graph end
#------------------------------------------------------------------------


    ###----------------------------------------------------------------------
    ### sub-function draw_solution start
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
            if  ( f_draw_monotone )
                mrp_ve_m = Morphing_monotonize( m_ve );
                p_ve_m = Morphing_extract_prm( mrp_ve_m );
                m_ve_m = Polygon_as_matrix( p_ve_m );
                plot!(plt, m_ve_m[1,:], m_ve_m[2,:],
                    linewidth=4,
                    label=:none, ticks=false, showaxis=false,
                    grid=:false, legend=false, framestyle=:none, lc=:yellow);
            end
        end
        f_debug  &&  println( "Drawing arrows..." );
    end
    ### sub-function draw_solution end
    #------------------------------------------------------------------------

    #------------------------------------------------------------------------
    function  draw_grid( plot )
        qlx = (ql[2:end-1])'
        plx = (pl[2:end-1])'

        if  f_draw_graph_only  &&  (! f_draw_grid )
            plot!( plt, [0; len_P], [qlx;qlx], lw=0.5, lc=:lightblue,
                   legend=false)
            plot!( plt,  [plx;plx], [0; len_Q], lw=0.5, lc=:lightblue,
                   legend=false)
        else
            plot!( plt, [0; len_P], [qlx;qlx], lw=2, lc=:black, legend=false)
            plot!( plt,  [plx;plx], [0; len_Q], lw=2, lc=:black, legend=false)
        end
    end
    ### sub-function draw_grid end
    #------------------------------------------------------------------------

    function  fz(x,y)
        p = Polygon_get_point_on( P, pl, x );
        q = Polygon_get_point_on( Q, ql, y );

        return Dist( p, q );
    end
    ### sub-function fz end
    #------------------------------------------------------------------------


    f_debug::Bool = false;

    f_debug  && println( "Getting ready to draw heatmap/graph/curves..." );

    len_P = total_length( P )
    len_Q = total_length( Q )

    pl = Polygon_prefix_lengths( P );
    ql = Polygon_prefix_lengths( Q );

    x_range = range(0, len_P, length = 200 )
    y_range = range(0, len_Q; length = 200 )


    f_debug && println( "Computing heat map..." );
    ### f_draw_graph_only
    if  f_draw_graph_only  && ( ! f_draw_grid )
        plt = plot( x_range, y_range, 0,
                    left_margin = 0 * Plots.mm,
                    bottom_margin=0*Plots.mm,
                    right_margin = 10.2 * Plots.mm,
                    ticks = false, showaxis = false, framestyle=:none,
                    dpi = 200 );
    else
        if   f_3d
            plt = plot( x_range, y_range, fz,
                       color = :haline,
                       left_margin = 0 * Plots.mm,
                       bottom_margin=0*Plots.mm,
                       right_margin = 10.02 * Plots.mm,
                        ticks = true, showaxis = true,
                        framestyle=:box,
#                        framestyle=:zerolines,
                       dpi = 300, st=:surface, camera=(30,50) );
        else
            plt = heatmap( x_range, y_range, fz,
                           color = :haline,
                           left_margin = 0 * Plots.mm,
                           bottom_margin=0*Plots.mm,
                           right_margin = 10.02 * Plots.mm,
                           ticks = false, showaxis = false, framestyle=:none,
                           dpi = 200 );
        end
    end

    if  ( length( title ) > 0 )
        title!( plt, title );
    end
    f_debug && println( "Heat map drawing done..." );



    cardi::Int64 = cardin(P) + cardin( Q );
    if  ( cardi < 2000 )
        if f_draw_c  ||  f_draw_ve  ||  f_draw_graph || f_draw_grid
            draw_grid( plt );
        end
    end
    if  ( cardi < 2000 )  &&  f_draw_graph
        f_debug  && println( "Drawing the graph..." );
        draw_graph( plt );
    else
        (! f_draw_grid)  &&  draw_solution( plt )
    end

    if  ( f_m_out_defined )
        p_c_diag = Morphing_extract_prm( m_out );
        m_c_diag = Polygon_as_matrix( p_c_diag );
        plot!(plt, m_c_diag[1,:], m_c_diag[2,:],
              linewidth=2, label=:none, ticks=false,
              showaxis=false, grid=:false,
              legend=false, framestyle=:none, lc=:red);
    end

    f_debug && println( "Saving heatmap/graph... ", filename_diagram );
    savefig( plt, filename_diagram );

    f_debug && println( "Outputing the curves..." );
#    output_polygons_to_file(  [P, Q], filename_curves, true );
    f_debug && println( "Created: " );
#    println( "   ", filename_curves );
    f_debug && println( "   ", filename_diagram );
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


function  html_open_w_file( filename::String, title::String )
    fl = open( filename, "w" )

    #println( "Writing file\n\n\n\n\n" );

    style = """
      <style>
      table {
         position: relative;
      }
      table::before,
      table::after {
         border: 1px solid #FFF;
         content: "";
         height: 100%;
         position: absolute;
         top: 0;
         width: 6px;
      }
      table::before {
         border-right: 0px;
         left: -6px;
      }
      table::after {
         border-left: 0px;
         right: -6px;
      }
      td {
         padding: 2px;
         text-align: center;
    }
    table, th, td {
    border: 1px outset black;
    }
    </style>
    """;

    write( fl, "<head>\n"
           * "<meta charset=\"UTF-8\">"
           * "<TITLE>" * title * "</TITLE>\n"
           * "<script type=\"text/x-mathjax-config\">\n"
           * "MathJax.Hub.Config({ tex2jax: "
           * "{inlineMath: [[\'\$\',\'\$\'], [\'\\(','\\)']]}"
           * "});\n"
           * "</script>\n"
           * "<script type=\"text/javascript\"\n"
           * "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js"
           * "config=TeX-AMS-MML_HTMLorMML\">\n"
           * "</script>\n"
           * "<meta charset=\"UTF-8\">\n"
           * style
           * "</head>" )
    write( fl, "<body>\n" );

    return  fl;
end


function  html_close( fl )
    println( fl, "<hr>\n" );
    dt=now();
    date_str = Dates.format(dt, "yyyy-mm-dd HH:MM:SS")
    write( fl, date_str );

    write( fl, "\n\n</body>\n" );
    close( fl );
end

function   html_write_video_file( fl, mvname::String,
    title::String )

    write( fl, "\n<h2>"*title*"</h2>\n\n" );

    write( fl, "<video controls autoplay " );
    write( fl, "   src=\""*mvname*"\" "
               * " type=\"video/mp4\" />\n" );
    write( fl, "</video>\n\n" );
end


function  draw_m_polygon( m::Morphing2F, curves_out )
    list = VecPolygon2F();

    push!( list, deepcopy( m.P ) );
    push!( list, deepcopy( m.Q ) );

    bb = BBox2F();
    BBox_bound( bb, m.P );
    BBox_bound( bb, m.Q );

    p = BBox_bottom_left( bb );
    for  poly  in list
        Polygon_translate!( poly, p );
    end

    output_polygons_to_file( list, curves_out * ".pdf", true, true );
    output_polygons_to_file( list, curves_out * ".png", false, true );
end




function  create_demo( title::String,
                       prefix::String,
                       poly_a, poly_b,
                       f_draw_c::Bool = false,
                       f_draw_ve::Bool = true,
                       note::String = "",
                       f_refinements::Bool = false )
    f_debug::Bool = false;

    println( "Creating: ", title );

    if  ! isdir( prefix )
        mkdir( prefix );
    end
    cardi = cardin( poly_a ) + cardin( poly_b );
    f_debug  &&  println( "Cardinality of both polygons: ", cardi );
    total_frames = min( 50 * (cardin( poly_a ) + cardin( poly_b )), 800 );

    filename_curves = prefix*"curves.pdf";
    f_debug && println( "Outputing the curves..." );
    output_polygons_to_file(  [poly_a, poly_b], filename_curves, true );

    options_svg = [ prefix*"curves.pdf", prefix*"curves.svg" ];
    output = read(pipeline( `pdf2svg $options_svg`, stderr="/tmp/errs.txt" ),
        String);

    println( "Created: ", filename_curves );

    local P, Q, m_d, m_d_r, m_ve_r, m_refinments;
    local m_d_dtw, m_SweepDist;

    m_SweepDist_vec = Vector{Morphing2F}();
    SweepDist_lb_vec = Vector{Float64}();

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
    f_debug && println( "A0: frechet_c_compute done..." );
    if  f_draw_ve
        f_debug && println( "A0: frechet_ve_r_compute about to be called..." );
        println( "Before m_ve_r..." );
        @time m_ve_r = frechet_ve_r_compute( poly_a, poly_b );
        println( "Iterations: ", m_ve_r.iters );
    end

    f_debug && println( "A1..." );

    f_SweepDist::Bool = false;
    if  ( cardi < 5000 )
        f_SweepDist = true;
        if  ( cardi < 100 )
            f_sampled_10 = true;
            P = Polygon_sample_uniformly( poly_a, 10*cardin( poly_a ) );
            Q = Polygon_sample_uniformly( poly_b, 10*cardin( poly_a ) );
        else
            P = poly_a;
            Q = poly_b;
        end
        f_debug && println( "A2..." );
        m_d = frechet_d_compute( P, Q );
        f_debug && println( "A3..." );
        m_d_dtw = DTW_d_compute( P, Q );
        f_debug && println( "A4..." );
        m_d_r = frechet_d_r_compute( P, Q );
        f_debug && println( "A5..." );
        f_computed_d = true;
    end

    if  f_SweepDist
        f_debug && println( "A6..." );
        m_SweepDist = SweepDist_compute( poly_a, poly_b );
        f_debug && println( "A7..." );
        m_SweepDist_r_m = SweepDist_compute_refine_mono( poly_a, poly_b );
        f_debug && println( "A8..." );

        SweepDist_compute_split( poly_a, poly_b, m_SweepDist_vec, 7, 2000 );

        for  i in eachindex( m_SweepDist_vec )
            mr = m_SweepDist_vec[ i ];
            mlb = SweepDist_lb_compute( mr.P, mr.Q )
            push!( SweepDist_lb_vec, mlb.sol_value );
#           exit( -1 );
        end
        Pa, Qa = Morphing_as_polygons( m_SweepDist_r_m );
        m_SweepDist_r_m_2 = SweepDist_compute( Pa, Qa );

        f_debug && println( "A9..." );
    end

    local m_refinements::Vector{Morphing2F} = Vector{Morphing2F}();
    #println( "BOGI!\n\n\n\n" );
    if   f_refinements
        m_refinements = Vector{Morphing2F}();
        frechet_mono_via_refinement_ext( poly_a, poly_b, m_refinements, true,
                                         1.000000001
                                      );
        f_debug && println( "m_refinements.len: ", length(m_refinements ) );
    end

    f_debug && println( "A10..." );

    #####################################################################
    # Creating movies/diagrams/etc
    #####################################################################

    output_polygons_to_file(  [poly_a, poly_b], prefix * "curves.png", false;
                          u_width=10.0);
    f_debug && println( "A10.a..." );

    f_debug && println( "A10.b..." );
    is_mkdir( prefix*"ortho/" );

    f_debug && println( "A11..." );

    if  f_draw_ve
        PU, QU = Morphing_as_polygons( m_ve_r );
        output_polygons_to_file(  [PU, QU], prefix*"ve_matching.pdf",
                                  true, true, true );

    end

    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.pdf",
        false, false, false
    );
    plot_curves_diagram( poly_a, poly_b, prefix*"diagram.png",
        false, false, false
    );
    plot_curves_diagram( poly_a, poly_b, prefix*"diagram_3d.png",
        false, false, false; f_3d=true
    );

    f_grid_only_drawn::Bool = false;
    if  ( cardi < 100 )
        plot_curves_diagram( poly_a, poly_b, prefix*"grid_only.pdf",
            false, false, false; f_draw_grid = true );
        plot_curves_diagram( poly_a, poly_b, prefix*"grid_only.png",
            false, false, false; f_draw_grid = true );
        f_grid_only_drawn = true;
    end

    f_debug && println( "A12..." );

    if   f_refinements
        dir =  prefix * "steps/";
        is_mkdir( dir );

        fl_s = html_open_w_file( dir * "index.html", "Refinement" );
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

            base_out = @sprintf( "curves_%06d", i );
            curves_out = dir * base_out;
            draw_m_polygon( mx, curves_out );

            title_frm = @sprintf( "Frame %02d   Monotonicity error: %s%%",
                              i, s )
            png_base = @sprintf( "%06d.png", i );
            png_out = dir * png_base;
            plot_curves_diagram( poly_a, poly_b,
                                 png_out,
                                 false, false, false, true,
                                 mx,
                                 title_frm
                                 );
            write( fl_s, "<hr>\n\n" );
            write( fl_s, "<img src=\"" * png_base * "\" />\n" )
            write( fl_s, "<hr>\n" )
            write( fl_s, "<img src=\"" * base_out * ".png" * "\" />\n" )
        end
        html_close( fl_s );

        println( "Generating gif..." );

        options = [ "-delay", "50", "-loop", "0", dir * "*.png",
                    prefix * "refinements.gif" ];
        outx = read(pipeline( `convert $options`, stderr="/tmp/errs.txt" ),
                    String );
    end

    f_graph_drawn::Bool = false;
    if  ( cardi < 500 )
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
        plot_curves_diagram( poly_a, poly_b, prefix*"ve_r_diagram_m.pdf",
            false, true, true;     f_draw_monotone = true );
        plot_curves_diagram( poly_a, poly_b, prefix*"ve_r_diagram_m.png",
            false, true, true; f_draw_monotone = true );
    end
    # println( "   ", prefix*"f_c_movie.mp4" );



    write_to_file( poly_a, prefix * "poly_a.txt" );
    write_to_file( poly_b, prefix * "poly_b.txt" );

    if  f_computed_d
        output_polygons_to_file(  [P, Q],
                                  prefix * "polygons_sampled.png", false,
                                  true  );
    end
    #########################################################
    # Html file...
    #########################################################
    fl = html_open_w_file( prefix*"index.html", prefix );

    write( fl, "<h1>", title, "</h1>\n" );


    #println( "Cardinality of both curves : ", cardi );



    # Image of curves

    write( fl, "<hr>\n\n" );
    write( fl, "<img src=\"curves.png\" />\n" )
    write( fl, "<hr>\n" )

    write( fl, "<h2><a href=\"SweepDist/\">SweepDist</a></h2>\n<hr>\n\n" );

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
              string( total_length( poly_a ) ) ]
    row_a_x = permutedims( row_a );

    row_b =  ["<a href=\"poly_b.txt\">Q</a>"
              string( cardin( poly_b ) )
              string( total_length( poly_b ) )
              ];
    row_b_x = permutedims( row_b );

    data = vcat( row_a_x, row_b_x  );

    style_tbl = Dict{String, String}();
    style_tbl[ "border" ] = "4px solid green";
    style_tbl[ "bordercolor" ] = "blue";
    pretty_table(fl, data;
                 header = (["Curves", "# Vertices", "Length"]),
                 table_style = style_tbl,
                 allow_html_in_cells = true,
        backend = Val(:html)
#        standalone=true
    )
    write( fl, "\n<hr>\n" )

    if  f_draw_ve
        B=fill("", (2,3) )
        B[1,1] = "Fréchet";
        B[1,2] = string( m_c.leash );
        B[1,3] = string( m_c.iters );
        B[2,1] = "VE Fréchet";
        B[2,2] = string( m_ve_r.leash ) ;
        B[2,3] = string( m_ve_r.iters );
        pretty_table(fl, B;
                     header = (["Distance",  "Value", "Iters"]),
                     allow_html_in_cells = true,
                     table_style = style_tbl,
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

    write( fl, "<h2>Free space diagram heatmap:</h2>\n\n" )
    write( fl, "<img src=\"diagram.png\" />\n" );


    if  ( f_graph_drawn )
        write( fl, "\n<hr>\n"
            *  "<a href=\"g_diagram.png\">graph + free space</a>\n"
            * "[<a href=\"g_diagram.pdf\">PDF</a>]\n"
            * " : "
            *  " <a href=\"graph_only.png\">graph</a>\n"
            * "[<a href=\"graph_only.pdf\">PDF</a>]\n<hr>\n"
        );
    end

    write( fl, "\n<hr>\n" )

    if  ( f_grid_only_drawn )
        write( fl, "<h3>With the grid</h3>\n\n" );

        write( fl, "<img src=\"grid_only.png\">\n" );

        write( fl, "\n<hr>\n" )
    end
    if  ( f_draw_ve  )
            write( fl, "<h2>VE-"*FrechetStr*" Retractable solution:</h2>" )
        write( fl, "<img src=\"ve_r_diagram.png\">\n" );
        write( fl, "<br><a href=\"ve_r_diagram_m.png\">Monotonized...</a>\n" );

    end

    write( fl, "<h2>" * FrechetStr * " cont+monotone solution:</h2>" )

    #        println( fl, "What the animation shows.<br>\n" );
    write( fl, "<img src=\"c_diagram.png\">\n" );
    write( fl, "<hr>\n" );


    if  ( f_computed_d )
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

    dir_SweepDist = prefix*"SweepDist/";
    local  fl_SweepDist;

    if  ( f_SweepDist )
        is_mkdir( dir_SweepDist )
        fl_SweepDist = html_open_w_file( dir_SweepDist * "index.html",
                                    "SweepDist distances/morhpings" );

        html_write_video_file( fl_SweepDist,  "SweepDist.mp4",
            "SweepDist (not necessarily monotone)" );

        html_write_video_file( fl_SweepDist, "SweepDist_r_m.mp4",
            "SweepDist monotonize via refinement" );


        for  i in eachindex( m_SweepDist_vec )
            str_num = @sprintf( "%06d", i )
            mvname = @sprintf( "%06d.mp4", i )
            #output_frechet_movie_mp4( m_SweepDist_vec[ i ],
            #                 dir_SweepDist * mvname, 400, true );
            html_write_video_file( fl_SweepDist, mvname,
                "Level "*str_num* " of SweepDist (splitting edges)<br>\n" );
        end

        println( fl_SweepDist, "<hr>\n" );
        A = fill("", 3+length( m_SweepDist_vec ), 6)
        for  i in eachindex( m_SweepDist_vec )
            m = m_SweepDist_vec[ i ];
            A[i,1] = string( i )
            A[i,2] = string( SweepDist_lb_vec[ i] );
            A[i,3] = string( Morphing_SweepDist_approx_price( m ) );
            A[i,4] = string( Morphing_SweepDist_price( m ) );
            A[i,5] = string( cardin( m.P ) );
            A[i,6] = string( cardin( m.Q ) );
            #=println( fl_SweepDist, i );
            println( fl_SweepDist, ": ", SweepDist_lb_vec[ i], "...",
                Morphing_SweepDist_approx_price( m ), " Exact: ",
                Morphing_SweepDist_price( m ) );
            println( fl_SweepDist, "<br>\n" );
            =#
            #            writeln( fl,SweepDist, "\n );
        end
        ind = length( m_SweepDist_vec ) + 1;
        A[ ind, 1 ] = "Sweep dist (orig): ";
        A[ ind, 4 ] =  string( Morphing_SweepDist_price( m_SweepDist ) );
        A[ ind, 5 ] = string( cardin( m_SweepDist.P ) );
        A[ ind, 6 ] = string( cardin( m_SweepDist.Q ) );

        A[ ind + 1, 1 ] = "Sweep dist r_mono: ";
        A[ ind + 1, 4 ] = string( Morphing_SweepDist_price( m_SweepDist_r_m ) );
        A[ ind + 1, 5 ] = string( cardin( m_SweepDist_r_m.P ) );
        A[ ind + 1, 6 ] = string( cardin( m_SweepDist_r_m.Q ) );

        A[ ind + 2, 1 ] = "Sweep dist r_mono_2: ";
        A[ ind + 2, 4 ] = string(
                 Morphing_SweepDist_price( m_SweepDist_r_m_2 ) );
        A[ ind + 2, 5 ] = string( cardin( m_SweepDist_r_m_2.P ) );
        A[ ind + 2, 6 ] = string( cardin( m_SweepDist_r_m_2.Q ) );

        println( fl_SweepDist, "<h2>Original SweepDist" );



        pretty_table(fl_SweepDist, A; header = (["Iteration", "Lower
                 bound", "Upper bound", "Better UB", "#P", "#Q" ]),
                 allow_html_in_cells = true, backend = Val(:html) );

        println( fl_SweepDist, "\n\n<hr>\n" );

        html_close( fl_SweepDist );
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
        write( fl, "<a href=\"steps/\">More info</a>\n\n" );
    end

    write( fl, "\n\n<h2>Animation: "*FrechetStr*" morphing as morphing</h2>" );
    write( fl, "\n\n <video controls autoplay " );
    write( fl, "   src=\"ortho/c.mp4\" type=\"video/mp4\" />\n" );
    write( fl, "</video>\n" );

    #=
    if  f_SweepDist
        println( fl, "<hr>" * "\n" );
        println( fl, "<h2>SweepDist</h2>\n" );
        write( fl, "\n\n <video controls autoplay " );
        write( fl, "   src=\"SweepDist.mp4\" type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );
        println( fl, "<h2>SweepDist refined monotone</h2>\n" );
        write( fl, "\n\n <video controls autoplay " );
        write( fl, "   src=\"SweepDist_r_m.mp4\" type=\"video/mp4\" />\n" );
        write( fl, "</video>\n" );
    end
    =#
    html_close( fl );

    println( "Generating all the movie files..." );
    if  f_draw_ve
        create_movie( poly_a, poly_b, total_frames,
                      prefix*"f_ve_r_movie.mp4", m_ve_r );
    end
    create_movie( poly_a, poly_b, total_frames, prefix*"f_c_movie.mp4", m_c );
    output_ortho_frechet_movie_mp4(  m_c, prefix*"ortho/c.mp4" );
    if  ( f_computed_d )
        output_frechet_movie_mp4( m_d_dtw, prefix*"d_dtw.mp4",
                                  400, true );
        output_frechet_movie_mp4( m_d, prefix*"discrete_frechet.mp4",
                                  400, true );
        output_frechet_movie( m_d, prefix*"discrete_frechet.pdf",
                              100, true );
        output_frechet_movie_mp4( m_d_r, prefix*"discrete_r_frechet.mp4",                                  400, true );
    end;
    if   f_SweepDist
        output_frechet_movie_mp4( m_SweepDist,
                                  dir_SweepDist * "SweepDist.mp4",
                                  400, true );
        output_frechet_movie_mp4( m_SweepDist_r_m,
                                  dir_SweepDist * "SweepDist_r_m.mp4",
                                  400, true );

        is_mkdir( dir_SweepDist )
        for  i in eachindex( m_SweepDist_vec )
            str_num = @sprintf( "%06d", i )
            mvname = @sprintf( "%06d.mp4", i )
            output_frechet_movie_mp4( m_SweepDist_vec[ i ],
                                      dir_SweepDist * mvname,
                                      400, true );
        end
    end

end


function  test_file( f_a )
    if  isfile( f_a )
        return
    end
    println( "\n\nERROR: Unable to open file!\n\n\t", f_a, "\n\n" );
    exit( -1 );
end


function  create_demo_files( title::String,
                             prefix::String,
                             f_a::String,
                             f_b::String,
                             f_draw_c::Bool = false,
                             f_draw_ve::Bool = true,
                             note::String = ""
                             )
    test_file( f_a );
    test_file( f_b );
    if ( ! isfile( f_a ) )
        return;
    end
    if ( ! isfile( f_b ) )
        return;
    end

    poly_a = polygon.read_file( f_a );
    poly_b = polygon.read_file( f_b );
    println( "#poly_a: ", cardin( poly_a ) );
    println( "#poly_b: ", cardin( poly_b ) );
    create_demo( title, prefix, poly_a, poly_b, f_draw_c, f_draw_ve, note );

end

gurl_base::String = "https://www.microsoft.com/en-us/research/publication/";
gurl_suffix::String = "geolife-gps-trajectory-dataset-user-guide/";

gurl::String = gurl_base * gurl_suffix;

gemb::String = "<a href=\"" * gurl * "\">GeoLife GPS Trajectories</a>";

function  gen_example_12()
    ( ! isfile( "data/geolife/5314.txt" ) )  &&  return;
    ( ! is_rebuild( "output/12" ) )  &&  return;

    create_demo_files( "Example of close curves (GPS tracks)",
                       "output/12/",
                       "data/geolife/5314.txt",
                       "data/geolife/5428.txt",
                       true, false,
                       "An example of two GPS tracks from " *
                           gemb * " that are close together. \n" *
                       "This is an example where the retractable Fréchet\n" *
                       " algorithm axplores only tiny fraction of the " *
                       " diagram, \n" *
                       "yielding a near linear running time in this case.\n"
                       );
end


function  gen_example_31()
    heart(t) = (-16*(sin(t))^3,
                13.0*cos(t)-5.0*cos(2*t) -2.0 *cos(3*t) - cos(4.0*t) )
    uheart(t) = (16*(sin(t))^3,
                 -(13.0*cos(t)-5.0*cos(2*t) -2.0 *cos(3*t) - cos(4.0*t) ) )

    poly_a = Polygon_fill( Polygon2F(), heart, -pi:0.02:pi );
    poly_b = Polygon_fill( Polygon2F(), uheart, 0:0.02:2.0*pi );
    create_demo( "Example 31:Inverse hearts",
                 "output/31/",
                 poly_a,poly_b,
                 true, true,
                 "Matching of two hearts",
                 true
                 );
end

function  gen_example_30()
    ( ! is_rebuild( "output/30" ) )  &&  return;

    poly_a,poly_b = example_30();
    create_demo( "Example 30: ZigZag does not math a Zag",
                 "output/30/",
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

function  gen_example_18()
    poly_a,poly_b = example_18();
    create_demo( "Example 18", "output/18/", poly_a, poly_b,
                 true, true,
        "Example where the minimum of the functions are linear",
                 true
                 );
end

function  gen_example_1()
    ( ! is_rebuild( "output/01" ) )  &&  return;

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
function  gen_example_32()
    poly_a,poly_b = example_32();
    create_demo( "Example 32", "output/32/",
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


function  rebuild_10()
    if  is_rebuild( "output/10" )
        println( "Example 10" );
        poly_a,poly_b = example_10( 3, 7);
        create_demo( "Example 10", "output/10/", poly_a,poly_b,
                     false, true, "", true );
    end
end

function  gen_example_2()
    ( ! is_rebuild( "output/02" ) )  &&  return;
    println( "Example 2" );
    poly_a,poly_b = example_2();
    create_demo( "Example 2", "output/02/", poly_a,poly_b );
end

function  gen_example_3()
    ( ! is_rebuild( "output/03" ) )  &&  return;
    println( "Example 3" );
    poly_a,poly_b = example_3();
    create_demo( "Example 3", "output/03/", poly_a,poly_b );
end

function  gen_example_4()
    ( ! is_rebuild( "output/04" ) )  &&  return;
    println( "Example 4" );
    poly_a,poly_b = example_4();
    create_demo( "Example 4", "output/04/", poly_a,poly_b );
end

function  gen_example_5()
    ( !is_rebuild( "output/05" ) )  &&  return;
    println( "Example 5" );
    poly_a,poly_b = example_5();
    create_demo( "Example 5", "output/05/", poly_a,poly_b,
                 false, true, "", true );
end

function  gen_example_6()
    ! is_rebuild( "output/06" )  &&  return;
    println( "Example 6" );
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

function  gen_example_33()
    ! is_rebuild( "output/33" )  &&  return;
    println( "Example 33" );
    poly_a,poly_b = example_10( 0, 7);
    create_demo( "Example 33", "output/33/", poly_a, poly_b,
                 false, true, "", true );
end


function  gen_example_34()
    ( ! is_rebuild( "output/34" ) )  &&  return;

    P, Q = example_34();

    create_demo( "Example 34", "output/34/",  P, Q,
                 false, true,
                 "Simplification removes unnecessary noise..."
                 );
end



function  gen_example_35()
    ( ! is_rebuild( "output/35" ) )  &&  return;

    P, Q = example_35();

    create_demo( "Example 35", "output/35/",  P, Q,
                 false, true,
                 "Zig zag zog..."
                 );
end



function  generate_examples()
    gen_example_1()
    gen_example_2();
    gen_example_3();
    gen_example_4();
    gen_example_5();
    gen_example_6();

    if  is_rebuild( "output/07" )
        println( "Example 7" );
        poly_a,poly_b = example_7();
        create_demo( "Example 7", "output/07/", poly_a,poly_b );
    end

    if  is_rebuild( "output/08" )
        println( "Example 8" );
        poly_a,poly_b = example_8_ext(3);
        create_demo( "Example 8 (3)", "output/08/", poly_a,poly_b );
    end

    if  is_rebuild( "output/09" )
        println( "Example 9" );
        poly_a,poly_b = example_9(  3, 4 );
        create_demo( "Example 9 (double zig-zag 3+4)", "output/09/",
                     poly_a,poly_b );
    end

    rebuild_10();

    if   is_rebuild( "output/11" )
        println( "Example 11" );
        # data/010/trajectory/20070910074631
        # data/010/trajectory/20070919122147
        create_demo_files( "Example of big curves (GPS tracks)",
                           "output/11/",
                           "data/geolife/8718.txt",
                           "data/geolife/1328.txt",
                           true, false );
    end

    gen_example_12();

    if  is_rebuild( "output/13" )
        println( "Example 13" );
        poly_a,poly_b = example_13();
        create_demo( "Example 13 (7 fixed)", "output/13/", poly_a,poly_b );
    end

    if  is_rebuild( "output/14" )
        println( "Example 14" );
        gen_example_14()
    end

    if  is_rebuild( "output/15" )
        println( "Example 15" );
        gen_example_15()
    end

    if  is_rebuild( "output/16" )
        println( "Example 16" );
        gen_example_16()
    end
    if  is_rebuild( "output/17" )
        println( "Example 17" );
        gen_example_17()
    end
    if  is_rebuild( "output/18" )
        println( "Example 18" );
        gen_example_18()
    end

    if   is_rebuild( "output/19" )
        println( "Example 19" );
        if  ( isfile( "data/birds/1787_1.plt" ) )
            create_demo_files( "Example of birds migration (GPS tracks)",
                               "output/19/",
                               "data/birds/1787_1.plt",
                               "data/birds/1787_2.plt",
                               true, true );
        end
    end

    if   is_rebuild( "output/20" )
        println( "Example 20" );
        if  ( isfile( "data/birds/1787_1.plt" ) )
            create_demo_files( "Example II of birds migration (GPS tracks)",
                "output/20/",
                "data/birds/2499_2.plt",
                "data/birds/2301_4.plt",
                true, true );
        end
    end

    if   is_rebuild( "output/21" )
        println( "Example 21" );
        if  ( isfile( "data/birds/2322_2.plt" ) )
            create_demo_files( "Example III of birds migration (GPS tracks)",
                               "output/21/",
                               "data/birds/2322_2.plt",
                               "data/birds/1793_4.plt",
                               true, true );
        end
    end


    if   is_rebuild( "output/22" )
        println( "Example 22" );
        create_demo_files( "Paper example 1 (birds: 1787_1 vs 1797_1",
                           "output/22/",
                           "data/birds/1787_1.plt",
                           "data/birds/1797_1.plt",
                           true, false );
    end
    if   is_rebuild( "output/23" )
        println( "Example 23" );
        create_demo_files( "Paper example 2 (birds: 2307_3 vs 2859_3",
                           "output/23/",
                           "data/birds/2307_3.plt",
                           "data/birds/2859_3.plt",
                           true, false );
    end
    if   is_rebuild( "output/24" )
        println( "Example 24" );
        create_demo_files( "Paper example 3 (birds: 2322_2 vs 1793_4",
                           "output/24/",
                           "data/birds/2322_2.plt",
                           "data/birds/1793_4.plt",
                           true, false );
    end


    if   is_rebuild( "output/25" )
        println( "Example 25" );
        # data/010/trajectory/20080928160000.plt
        # data/010/trajectory/20081219114010.plt
        create_demo_files(
            "Paper example 4: GeoLife 20080928160000 / 20081219114010",
            "output/25/",
            "data/geolife/1367_109.txt",
            "data/geolife/1403.txt",
            true, false );
    end
    if   is_rebuild( "output/26" )
        println( "Example 26" );
        # data/041/trajectory/20090708221430.plt
        # data/041/trajectory/20090712044248.plt
        create_demo_files(
            "Paper example 5: GeoLife 20090708221430 / 20090712044248",
            "output/26/",
            "data/geolife/5578.txt",
            "data/geolife/5585.txt",
            true, false );
    end
    if   is_rebuild( "output/27" )
        println( "Example 27" );
        create_demo_files(
            "Paper example 6: Pigeons RH887_1 / RH887_11",
            "output/27/",
            "data/pigeons/RH887_1.plt",
            "data/pigeons/RH887_11.plt",
            true, true );
    end
    if   is_rebuild( "output/28" )
        println( "Example 28" );
        create_demo_files(
            "Paper example 7: Pigeons C369_5 / C873_6",
            "output/28/",
            "data/pigeons/C369_5.plt",
            "data/pigeons/C873_6.plt",
            true, false );
    end
    if   is_rebuild( "output/29" )
        println( "Example 29" );
        create_demo_files(
            "Paper example 8: Pigeons C360_10 / C480_9",
            "output/29/",
            "data/pigeons/C360_10.plt",
            "data/pigeons/C480_9.plt",
            true, true );
    end


    if   is_rebuild( "output/30" )
        println( "Example 30" );
        gen_example_30()
    end

    if   is_rebuild( "output/31" )
        println( "Example 31" );
        gen_example_31()
    end

    if   is_rebuild( "output/32" )
        println( "Example 32" );
        gen_example_32()
    end

    gen_example_33();
    gen_example_34();
    gen_example_35();


end

function  reportIt( m )
    mlb = SweepDist_lb_compute( m.P, m.Q );
    println( mlb.sol_value, "...", Morphing_SweepDist_price( m ) );
end

function  SW_test()
    P, Q = example_3();

    for i in 1:80
        m = SweepDist_compute_refine_mono( P, Q );
        reportIt( m );

        i,j, mx = Morphing_get_max_edges_err( m );
        println( "Max edge fidelity: ", mx, " (", i, ", ", j, ")" );
#        PA = Polygon_split_single_edge( m.P, i+1 );
#        QA = Polygon_split_single_edge( m.Q, j+1 );

        l_p = Polygon_edge_length( m.P, i );
        l_q = Polygon_edge_length( m.Q, j );

        if  l_p > l_q
            P = Polygon_split_single_edge( m.P, i );
            Q = m.Q;
        else
            P = m.P;
            Q = Polygon_split_single_edge( m.Q, j );
        end
    end
end


####################################################################

test_command( "HandBrakeCLI", "ffmpeg", "pdf2svg", "convert" );

if  ! isdir( "output" );
    mkdir( "output" );
end

num_args = length( ARGS );

if   num_args == 1  &&  (ARGS[1]=="SD" )
    SW_test();
end

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

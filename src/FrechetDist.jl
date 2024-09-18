# Originally contributed by S. Har-Peled
# under MIT License

module FrechetDist

using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using DataStructures
using Printf


#include( "point.jl" );
#include( "polygon.jl" );

include( "cg.jl" );

using .cg
using .cg.point
using .cg.segment
using .cg.polygon


#using .Point;
#using .Polygon;
#using .Point
#using .Polygon
#using .cg
#using  cg;  # A submodule?

include( "morphing.jl" );
include( "frechet_discrete.jl" );
include( "frechet.jl" );
include( "SweepDist.jl" );
include( "utils.jl" );
include( "polygon_hierarchy.jl" );
include( "palette.jl" )
include( "fever.jl" )

##########################################################################




#-- Types -----------------------------------------------------
export  Morphing, FPointType, EventPoint#, PointType
export  Morphing2F

#export  FrechetCExtraInfo;

# Morphing stuff
export  Morphing_monotonize, Morphing_empty, Morphing_verify_valid
export  Morphing_combine, Morphing_extract_prm
export  Morphing_as_polygons
export  Morphing_extract_vertex_radii, Morphing_extract_offsets
export  Morphing_monotone_leash

export  frechet_refinement;
export  frechet_mono_via_refinement_delta


export  Morphing_SweepDist_price
export  Morphing_SweepDist_approx_price

export  Morphing_get_max_edges_err;

export  Morphing_sample_uniformly


#####################################################################
export  frechet_d_compute
export  frechet_d_r_compute
export  frechet_d_r_compute_sample

################
export  frechet_ve_r_compute
export  frechet_width_approx
export  frechet_offsets

export  frechet_ve_r_compute_range

export  FEVER_compute_range


export  frechet_mono_via_refinement, frechet_mono_via_refinement_ext

export  frechet_ve_r_compute_ext

export  frechet_c_mono_approx_subcurve;

export  frechet_dist_upper_bound

export  frechet_simplify_to_width
export  frechet_simplify_to_cardin
export  frechet_simplify_w_exp

export  frechet_c_approx

export  frechet_c_compute

export  DTW_d_compute;

export  SweepDist_compute;
export  SweepDist_lb_compute;
export  SweepDist_compute_refine_mono;
export  SweepDist_compute_split;

#########################################################################
# Palette: Is an array of floats that encodes approximation info about
# a polygon. Using such a precomputed palette, one can quickly extract
# a simplification.

export  frechet_palette
export  frechet_palette_approx
export  frechet_palette_level

#########################################################################
# Polygon hierarchy...
export   PolygonHierarchy

export  ph_push!
export  ph_print
export  ph_approx
export  ph_compute_hierarchy
export  ph_init




end

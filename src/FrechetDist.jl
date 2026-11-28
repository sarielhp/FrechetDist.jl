# Originally contributed by S. Har-Peled
# under MIT License

module FrechetDist

using BasicCompGeometry: BasicCompGeometry, AbsPntSeq, AnyDim, GPntSeq, PntSeq, PntSeq2F
using BasicCompGeometry: Point, Point2F
using BasicCompGeometry: cardin, dim, geom_length, point, point_max, prefix_lengths
using BasicCompGeometry: push_distinct!, split_edges
using BasicCompGeometry.ReadWrite: ReadWrite, read_plt
using DataStructures: DataStructures, BinaryHeap, BinaryMinHeap, PriorityQueue, SortedSet
using DataStructures: SetSemiToken, advance, beforestartsemitoken, deref, pastendsemitoken
using DataStructures: push_return_semitoken!, regress
using Random: Random, AbstractRNG
using LinearAlgebra: LinearAlgebra, lq, norm, qr, rank
using Parameters: Parameters, @with_kw
using StaticArrays: StaticArrays, MVector, StaticVector

# Extended for the event ordering of FEVER (fever.jl).
import Base.Order: lt

include("fp_geometry.jl")
include("fp_geometry_vec.jl")
include("curve_utils.jl")
include("morphing.jl")
include("morphing_prm.jl")
include("frechet_discrete.jl")
include("frechet.jl")
include("frechet_width.jl")
include("frechet_refine.jl")
include("frechet_c.jl")
include("frechet_decide.jl")
include("frechet_decide_blocks.jl")
include("frechet_decide_lazy_tree.jl")
include("frechet_decide_sausage.jl")
include("frechet_exact.jl")
include("frechet_bisect.jl")
include("frechet_simplify.jl")
include("SweepDist.jl")
include("utils.jl")
include("palette.jl")
include("polygon_hierarchy.jl")
include("fever.jl")
include("validate.jl")

##########################################################################




#-- Types -----------------------------------------------------
export Morphing, FPointType, EventPoint
export Morphing2F


# Morphing stuff
export Morphing_monotonize, Morphing_empty, Morphing_verify_valid
export Morphing_combine, Morphing_extract_prm
export Morphing_as_polygons
export Morphing_extract_vertex_radii, Morphing_extract_offsets
export Morphing_monotone_leash
export morphing_profile
export Morphing_as_function_w_times
export Morphing_SweepDist_price
export Morphing_SweepDist_approx_price
export Morphing_get_max_edges_err
export Morphing_sample_uniformly

export frechet_refinement
export frechet_mono_via_refinement_delta


export polygons_get_loc_at_time

#####################################################################
export frechet_d_compute
export frechet_d_r_compute
export frechet_d_r_compute_sample

################
export frechet_ve_r_compute
export frechet_width_approx
export frechet_offsets

export frechet_ve_r_compute_range

export FEVER_compute_range, FEVER_compute


export frechet_mono_via_refinement, frechet_mono_via_refinement_ext!

export frechet_ve_r_compute_ext

export frechet_c_mono_approx_subcurve

export frechet_dist_upper_bound

export frechet_simplify_to_width
export frechet_simplify_to_cardin
export frechet_simplify_w_exp

export frechet_c_approx

export frechet_c_compute
export frechet_c_decide_exact
export frechet_c_compute_exact
export frechet_c_compute_bisect

export DTW_d_compute

export SweepDist_compute
export SweepDist_lb_compute
export SweepDist_compute_refine_mono
export SweepDist_compute_split!

#########################################################################
# Palette: an array of floats that encodes approximation info about
# a polygon. Using such a precomputed palette, one can quickly extract
# a simplification.

export frechet_palette
export frechet_palette_approx
export frechet_palette_level

#########################################################################
# Polygon hierarchy...
export PolygonHierarchy

export ph_push!
export ph_print
export ph_approx
export ph_compute_hierarchy
export ph_init


export simplify_morphing_sensitive

end

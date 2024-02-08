module FrechetDist

using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using DataStructures
using Printf


include( "cg.jl" );
using .cg
#using .cg
#using  cg;  # A submodule?

include( "frechet.jl" );

##########################################################################




#-- Types -----------------------------------------------------
export Morphing, FPointType, EventPoint, PointType

# Morphing stuff
export  Morphing_monotonize, Morphing_empty, Morphing_verify_valid
export  Morphing_combine, Morphing_extract_prm
export  Morphing_as_polygons
export  Morphing_extract_vertex_radii, Morphing_extract_offsets

################
export  frechet_ve_r_compute
export  frechet_d_compute
export  frechet_d_r_compute
export  frechet_width_approx
export  frechet_offsets

export  frechet_mono_via_refinement

export  frechet_ve_r_compute_ext

export  frechet_dist_upper_bound

export  frechet_c_approx
export  frechet_c_compute


# Write your package code here.

end

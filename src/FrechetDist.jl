# Originally contributed by S. Har-Peled
# under MIT License

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

include( "morphing.jl" );
include( "frechet_discrete.jl" );
include( "frechet.jl" );

##########################################################################




#-- Types -----------------------------------------------------
export  Morphing, FPointType, EventPoint, PointType
export  Morphing2F

#export  FrechetCExtraInfo;

# Morphing stuff
export  Morphing_monotonize, Morphing_empty, Morphing_verify_valid
export  Morphing_combine, Morphing_extract_prm
export  Morphing_as_polygons
export  Morphing_extract_vertex_radii, Morphing_extract_offsets
export  Morphing_adtw_price


#####################################################################
export  frechet_d_compute
export  frechet_d_r_compute
export  frechet_d_r_compute_sample

################
export  frechet_ve_r_compute
export  frechet_width_approx
export  frechet_offsets

export  frechet_mono_via_refinement, frechet_mono_via_refinement_ext

export  frechet_ve_r_compute_ext

export  frechet_dist_upper_bound

export  frechet_c_approx
export  frechet_c_compute

export  DTW_d_compute;

export  ADTW_compute;
export  ADTW_compute_refine_mono;
export  ADTW_compute_split;

end

# Hierarchy of simplifications of a polygon, at decreasing widths. The
# hierarchy itself lives in BasicCompGeometry.PntSeqHierarchy; here we
# only define the Frechet measure controlling it.

"""
    frechet_measure

The measure of a polygon hierarchy under the Frechet distance: the
width of a subcurve is approximated greedily, simplification is the
greedy Frechet simplification, and the quality of a simplification is
certified by computing a monotone matching to it. Palettes use the
midpoint splitting rule.
"""
const frechet_measure = PH.HierarchyMeasure(
    frechet_width_approx,
    frechet_simplify_w_exp,
    (P, R, indices) -> frechet_c_mono_approx_subcurve(P, R, indices)[1].leash,
    PH.MidpointSplit(),
);

const PolygonHierarchy = PH.SimplificationHierarchy

const ph_push! = PH.hierarchy_push!
const ph_print = PH.hierarchy_print
const ph_approx = PH.hierarchy_approx
const ph_compute_hierarchy = PH.hierarchy_build!

function ph_init(P::PntSeq{D,T}, ph_lock) where {D,T}
    check_curve(P, "ph_init", "the curve")
    ph = PH.hierarchy_init(P, frechet_measure; lock = ph_lock)
    ph.len = fp_geom_length(P)

    return ph
end

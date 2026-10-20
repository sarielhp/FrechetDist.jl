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
const ph_compute_hierarchy = PH.hierarchy_build!

"""
    ph_approx( ph, w, n_max = 25, resolution = 4.0 )

An approximation of the hierarchy's curve with Frechet error at most about `w`,
returned as `(curve, width)`. The hierarchy caches the simplifications it computes, so
the result can depend on the calls made before. The returned curve is a copy: changing
it does not affect the hierarchy.
"""
function ph_approx(ph::PolygonHierarchy, w::Float64, args...)
    Q, width = PH.hierarchy_approx(ph, w, args...)
    return typeof(Q)(copy(Q.pnts)), width
end

"""
    ph_init( P; lock = ReentrantLock() )

A polygon hierarchy of `P`. The lock guards the hierarchy when it is shared between
threads. (`ph_init( P, lock )` is accepted as well.)
"""
function ph_init(P::PntSeq{D,T}; lock::ReentrantLock = ReentrantLock()) where {D,T}
    check_curve(P, "ph_init", "the curve")
    ph = PH.hierarchy_init(P, frechet_measure; lock)
    ph.len = fp_geom_length(P)

    return ph
end

ph_init(P::PntSeq, ph_lock::ReentrantLock) = ph_init(P; lock = ph_lock)

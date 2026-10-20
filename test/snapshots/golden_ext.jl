# Supplemental snapshot for functions the parent golden.jl does NOT cover (found by audit of golden.jl):
# FEVER_compute(_range), frechet_ve_r_compute_range, frechet_d_r_compute, SweepDist_lb_compute,
# SweepDist_compute_refine_mono, frechet_c_approx (several ratios), frechet_mono_via_refinement(_ext/_delta),
# extract_refined_polygon(delta), Morphing_combine/prm round-trip, frechet_c_compute(f_accept_approx=false).
# Usage: julia --project=<proj> golden_ext.jl out.jls
using Random, Serialization
include(joinpath(@__DIR__, "..", "agreement.jl"))
import BasicCompGeometry: cardin, Points
const out = Dict{String,Any}()
mat(P) = [Points(P)[i][d] for d = 1:length(Points(P)[1]), i = 1:cardin(P)]
ev(e) = (Vector(e.p), e.i, Int(e.type), e.t, e.penalty)
conv(m::FD.Morphing) = (
    m.leash,
    m.iters,
    m.ratio,
    m.lower_bound,
    m.monotone_err,
    mat(m.P),
    mat(m.Q),
    ev.(m.pes),
    ev.(m.qes),
)
conv(t::Tuple) = map(conv, t)
conv(x::PntSeq) = mat(x)
conv(x) = x
function rec(k, f)
    out[k] = try
        redirect_stdout(devnull) do
            conv(f())
        end
    catch e
        "ERR " * first(sprint(showerror, e), 60)
    end
end
rng = MersenneTwister(2024)
for (fam, gen) in pair_families(), c = 1:4
    Pv, Qv = gen(rng)
    P = mkseq(Pv)
    Q = mkseq(Qv)
    (cardin(P) < 2 || cardin(Q) < 2) && continue
    id = "$fam/$c"
    ve = FD.frechet_ve_r_compute(P, Q)
    rec("$id/fever", () -> FD.FEVER_compute(P, Q))
    for f in (0.5, 1.0, 2.0)
        rec("$id/ve_range/$f", () -> FD.frechet_ve_r_compute_range(P, Q, f * ve.leash))
        rec("$id/fever_range/$f", () -> FD.FEVER_compute_range(P, Q, f * ve.leash))
    end
    rec("$id/d_r", () -> FD.frechet_d_r_compute(P, Q))
    rec("$id/sweep_lb", () -> FD.SweepDist_lb_compute(P, Q))
    rec("$id/sweep_refine", () -> FD.SweepDist_compute_refine_mono(P, Q))
    for a in (1.5, 3.0)
        rec("$id/c_approx/$a", () -> FD.frechet_c_approx(P, Q, a))
    end
    rec("$id/mono_ref", () -> FD.frechet_mono_via_refinement(P, Q)[1])
    rec(
        "$id/mono_ref_ext",
        () -> (
            o = Vector{FD.Morphing2F}();
            r = FD.frechet_mono_via_refinement_ext(P, Q, o, true, 1.001);
            (r[1], r[2], r[3], r[4], Tuple(o))
        ),
    )
    rec("$id/mono_delta", () -> FD.frechet_mono_via_refinement_delta(P, Q, 0.01, true))
    rec(
        "$id/refine_delta",
        () -> (
            FD.extract_refined_polygon(P, ve.pes, 2, 0.05),
            FD.extract_refined_polygon(Q, ve.qes, 1),
        ),
    )
    rec("$id/c_compute_noacc", () -> FD.frechet_c_compute(P, Q, false))
    rec(
        "$id/prm_roundtrip",
        () -> (
            mo = FD.Morphing_monotonize(ve);
            prm = FD.Morphing_extract_prm(mo);
            (
                FD.Morphing_from_prm(prm, P, Q),
                FD.Morphing_combine(mo, FD.Morphing_init(Q, P, mo.qes, mo.pes)),
            )
        ),
    )
    rec("$id/offsets", () -> FD.Morphing_extract_offsets(ve))
end
serialize(ARGS[1], out);
println("wrote ", length(out), " entries; errors: ", count(v -> v isa String, values(out)))

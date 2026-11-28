# Decider benchmark for FrechetDist on the decider queries of Bringmann, Künnemann and
# Nusser (the counterpart of decider_bench.cpp).
#
# Usage: julia --project=. examples/cpp_compare/fd_bench.jl <mode> <curve_dir> \
#            <answers_out> <runs> <query_file>...
#
# <mode> is `pid`, the decider of examples/run_queries_test.jl (frechet_decider_PID,
# copied below unchanged apart from the debug output), or `exact_lazy`, the library's
# exact decision procedure frechet_c_decide_exact on curves prepared when read with the
# default DecideCurve (checked, with their width and bounding box, as the C++ code
# computes the bounding box and prefix lengths of a curve when it reads it; the prefix
# lengths and the box hierarchy are built when first needed, inside the timed decisions:
# the curves are prepared again, untimed, before each timed run, so that each run builds
# what it needs). `exact` prepares them with DecideCurve(P, BoxTree), the hierarchy
# built when read (untimed); `exact_plain` calls frechet_c_decide_exact on the curves
# themselves. `pid:<eps>` is the PID decider with tolerance
# <eps>: it may answer either way when the distance is within a factor 1 +- eps of the
# query distance (it answers "yes" once its upper bound is below (1 + eps) r and "no"
# once its lower bound exceeds (1 - eps) r); `pid` is `pid:0`. `approx:<eps>` decides
# with frechet_c_approx(P, Q, 1 + eps): "yes" when its leash is below r.
#
# As in decider_bench.cpp, the curves of a query file are read first (with their
# simplification hierarchies for `pid`, built as in the serial mode of
# run_queries_test.jl, and prepared for `exact`), and only the decision calls are
# timed. Each file is run once untimed (warm-up, which also compiles), then timed
# <runs> times; the minimum is reported. Prints
# "<file> <queries> <yes> <min_ms> <prep_ms>" per query file and writes
# "<curve1> <curve2> <distance> <0|1>" per query (1 when the distance is below the
# query distance).
using FrechetDist, BasicCompGeometry, Printf
import FrechetDist: read_file
import BasicCompGeometry: cardin

const Lazy = FrechetDist.LazyBoxTree{Vector{Point2F}}

struct Curves
    polys::Vector{PntSeq2F}
    prepared::Vector{FrechetDist.DecideCurve{PntSeq2F,FrechetDist.BoxTree}}
    lazy::Vector{FrechetDist.DecideCurve{PntSeq2F,Lazy}}
    PHA::Vector{PolygonHierarchy}
    index::Dict{String,Int}
end

function load_curve!(C::Curves, dir, name, mode)
    f_hier = startswith(mode, "pid")
    haskey(C.index, name) && return C.index[name]
    P = read_file(joinpath(dir, name))
    push!(C.polys, P)
    mode == "exact" && push!(C.prepared, FrechetDist.DecideCurve(P, FrechetDist.BoxTree))
    mode == "exact_lazy" &&
        push!(C.lazy, FrechetDist.DecideCurve(P, FrechetDist.LazyBoxTree))
    if f_hier
        ph = ph_init(P, ReentrantLock())
        ph_compute_hierarchy(ph, 40.0, min(15, round(Int64, cardin(P) / 4)))
        push!(C.PHA, ph)
    end
    return C.index[name] = length(C.polys)
end

# frechet_decider_PID of examples/run_queries_test.jl: sign(d_F(P, Q) - r).
function frechet_decider_PID(
    PID::Curves,
    i::Int64,
    j::Int64,
    r::Float64,
    eps::Float64 = 0.0,
)::Int64
    r_lo = r * (1.0 - eps)  # answer "no" once the lower bound exceeds r_lo
    r_hi = r * (1.0 + eps)  # answer "yes" once the upper bound is below r_hi
    PA = P_orig = PID.polys[i]
    QA = Q_orig = PID.polys[j]

    l_b = dist(first(PA), first(QA))
    (l_b > r_lo) && return 1

    l_b = max(l_b, dist(last(PA), last(QA)))
    (l_b > r_lo) && return 1

    P_ph = PID.PHA[i]
    Q_ph = PID.PHA[j]

    wP = P_ph.widths[1]
    wQ = Q_ph.widths[1]

    ub = l_b + wP + wQ
    (ub < r_hi) && return -1

    delta_naive = (ub - l_b) / 4.0
    delta = min(abs(r - l_b), abs(r - ub), (wP + wQ)/4.0, delta_naive)
    if (delta < delta_naive / 10.0)
        delta = delta_naive / 10.0
    end
    if (delta < eps * r / 4.0)  # no finer than the tolerance needs
        delta = eps * r / 4.0
    end

    P_limit::Int64 = round(Int64, 0.82 * cardin(P_orig))
    Q_limit::Int64 = round(Int64, 0.82 * cardin(Q_orig))

    f_monotone::Bool = false
    f_orig::Bool = false

    for iters::Int64 in 1:20
        w_trg = delta / max(2.0, 1.0 + iters)
        if f_orig
            PA = P_orig
            QA = Q_orig
            wQ = wP = 0.0
        else
            PA, wP = ph_approx(P_ph, w_trg)
            QA, wQ = ph_approx(Q_ph, w_trg)
            if ((cardin(PA) > P_limit) || (cardin(QA) > Q_limit))
                QA = Q_orig
                PA = P_orig
                wP = wQ = 0.0
                f_orig = true
            end
        end

        if f_monotone
            m, PA_A, QA_A = frechet_mono_via_refinement_delta(PA, QA, delta, false)
            l_min = m.lower_bound
            l_max = m.leash
        else
            l_min, l_max = FEVER_compute_range(PA, QA, ub)
            if ((iters > 0) && (l_min < r < l_max) && ((l_max - l_min) > 1.1*delta))
                f_monotone = true
                delta = delta / 6.0
                continue
            end
        end
        l_b = max(l_b, l_min - wP - wQ)
        (l_b > r_lo) && return +1

        ub = min(ub, l_max + wP + wQ)
        (ub < r_hi) && return -1

        old_delta = delta
        delta = min(abs(l_max - r), abs(l_min - r), ub - r, r - l_b, delta / 2.0)
        if (delta < old_delta / 10.0)
            delta = old_delta / 10.0
        end
        if (delta < eps * r / 4.0)
            delta = eps * r / 4.0
        end
    end

    m = frechet_c_compute(P_orig, Q_orig)
    return round(Int64, sign(m.leash - r))
end

mode_eps(mode) = occursin(':', mode) ? parse(Float64, split(mode, ":")[2]) : 0.0

# The approximate decider: frechet_c_approx gives a leash d' with d <= d' <= (1 + eps) d,
# and the answer is whether d' < r. It is right unless (1 - eps) r < d <= r.
function decide_approx(P, Q, r, eps)
    max(dist(first(P), first(Q)), dist(last(P), last(Q))) > r && return 0
    return frechet_c_approx(P, Q, 1.0 + eps).leash < r ? 1 : 0
end

# The same answer, starting from a coarse approximation: a leash d' with ratio a gives
# d' / a <= d <= d', which decides whenever r is outside [d' / a, d']; otherwise the
# ratio shrinks towards what the gap between r and the bounds needs, down to 1 + eps.
function decide_approx_adaptive(P, Q, r, eps)
    max(dist(first(P), first(Q)), dist(last(P), last(Q))) > r && return 0
    ratio = 4.0
    while true
        m = frechet_c_approx(P, Q, ratio)
        m.leash < r && return 1
        lb = m.leash / m.ratio
        lb > r && return 0
        ratio <= 1.0 + eps && return 0
        need = max(m.leash / r, r / lb)  # a ratio below this separates r from a bound
        ratio = max(1.0 + eps, min(sqrt(ratio), 1.0 + (need - 1.0) / 2.0))
    end
end

# 1 when d_F(P, Q) < r, 0 when it is larger; -1 when the decider reports equality.
function decide(mode, C::Curves, i, j, r)
    if startswith(mode, "pid")
        s = frechet_decider_PID(C, i, j, r, mode_eps(mode))
        return s < 0 ? 1 : (s > 0 ? 0 : -1)
    end
    startswith(mode, "approx:") &&
        return decide_approx(C.polys[i], C.polys[j], r, mode_eps(mode))
    startswith(mode, "approx_adaptive:") &&
        return decide_approx_adaptive(C.polys[i], C.polys[j], r, mode_eps(mode))
    mode == "exact_plain" &&
        return frechet_c_decide_exact(C.polys[i], C.polys[j], r) ? 1 : 0
    mode == "exact_lazy" && return frechet_c_decide_exact(C.lazy[i], C.lazy[j], r) ? 1 : 0
    return frechet_c_decide_exact(C.prepared[i], C.prepared[j], r) ? 1 : 0
end

function bench_file(mode, dir, file, runs, answers)
    C = Curves(
        PntSeq2F[],
        FrechetDist.DecideCurve{PntSeq2F,FrechetDist.BoxTree}[],
        FrechetDist.DecideCurve{PntSeq2F,Lazy}[],
        PolygonHierarchy[],
        Dict{String,Int}(),
    )
    queries = Tuple{String,String,Int,Int,Float64}[]
    prep = @elapsed for line in eachline(file)
        f = split(line)
        length(f) == 3 || continue
        i = load_curve!(C, dir, f[1], mode)
        j = load_curve!(C, dir, f[2], mode)
        push!(queries, (f[1], f[2], i, j, parse(Float64, f[3])))
    end
    result = zeros(Int, length(queries))
    # Warm-up: one untimed pass, so that everything is compiled before the timed runs.
    for (k, (_, _, i, j, r)) in enumerate(queries)
        result[k] = decide(mode, C, i, j, r)
    end
    best = Inf
    for _ in 1:runs
        mode == "exact_lazy" &&
            map!(P -> FrechetDist.DecideCurve(P, FrechetDist.LazyBoxTree), C.lazy, C.polys)
        t = @elapsed for (k, (_, _, i, j, r)) in enumerate(queries)
            result[k] = decide(mode, C, i, j, r)
        end
        best = min(best, t)
    end
    for (k, (a, b, _, _, r)) in enumerate(queries)
        println(answers, a, " ", b, " ", repr(r), " ", result[k])
    end
    @printf(
        "%s %d %d %.3f %.3f\n",
        file,
        length(queries),
        count(==(1), result),
        1000best,
        1000prep
    )
    flush(stdout)
end

function main(args)
    length(args) >= 5 || error(
        "usage: fd_bench.jl <pid|pid:<eps>|approx:<eps>|exact|exact_plain> <curve_dir> <answers_out> <runs> <query_file>...",
    )
    mode, dir, out, runs = args[1], args[2], args[3], parse(Int, args[4])
    (
        mode == "exact" ||
        mode == "exact_plain" ||
        mode == "exact_lazy" ||
        startswith(mode, "pid") ||
        startswith(mode, "approx")
    ) || error("mode must be pid, pid:<eps>, approx:<eps> or exact")
    open(out, "w") do answers
        for file in args[5:end]
            bench_file(mode, dir, file, runs, answers)
        end
    end
end

main(ARGS)

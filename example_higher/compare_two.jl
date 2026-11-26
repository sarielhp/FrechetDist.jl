#!/usr/bin/env julia
# Compares two motion-capture trials, as curves of joint angles in ℝ^46, under the
# distances the package computes, and checks the order they must satisfy.
#
# Usage (from the repository root, after `ruby example_higher/download.rb`):
#   julia example_higher/compare_two.jl [ID1 ID2 [STRIDE]]
# IDs are trial ids from trials.txt (default 35_01 35_17: a walk and a run/jog of the
# same subject). STRIDE keeps every STRIDE-th frame of the 120 Hz recording (default 1).

using QuickEnv # local
using FrechetDist
using Printf

include(joinpath(@__DIR__, "mocap.jl"))
using .Mocap

const DISTANCES = [
    ("Fréchet (exact)", (P, Q) -> frechet_c_compute(P, Q).leash),
    ("Fréchet, 1.01-approximation", (P, Q) -> frechet_c_approx(P, Q, 1.01).leash),
    ("VE-retractable Fréchet", (P, Q) -> frechet_ve_r_compute(P, Q).leash),
    ("discrete Fréchet", (P, Q) -> frechet_d_compute(P, Q).leash),
    ("DTW (sum of distances)", (P, Q) -> DTW_d_compute(P, Q).sol_value),
]

function describe(t::Trial, n::Int)
    @printf(
        "  %s  %-5s  %4d frames (%d dropped), %4d vertices used\n",
        t.id,
        t.activity,
        length(t) + t.dropped,
        t.dropped,
        n
    )
end

# Runs every distance twice and times the second run, so compilation is not counted.
function measure(P, Q)
    return map(DISTANCES) do (name, f)
        f(P, Q)
        seconds = @elapsed value = f(P, Q)
        @printf("  %-30s %12.4f   %8.4f s\n", name, value, seconds)
        name => value
    end
end

function compare(a::Trial, b::Trial, ::Val{D}, stride::Int) where {D}
    P = trial_curve(a, Val(D); stride)
    Q = trial_curve(b, Val(D); stride)
    println("Curves in ℝ^$D (joint angles in degrees):")
    describe(a, length(P))
    describe(b, length(Q))
    println()
    values = Dict(measure(P, Q))
    exact = values["Fréchet (exact)"]
    checks = [
        "VE-retractable ≤ exact" =>
            values["VE-retractable Fréchet"] <= exact * (1 + 1e-9),
        "exact ≤ 1.01-approximation ≤ 1.01·exact" =>
            exact <= values["Fréchet, 1.01-approximation"] * (1 + 1e-9) <= 1.01 * exact,
        "exact ≤ discrete" => exact <= values["discrete Fréchet"] * (1 + 1e-9),
    ]
    println()
    for (what, ok) in checks
        println("  ", ok ? "ok    " : "FAILED", " ", what)
    end
    return all(last, checks)
end

function main(args)
    ids = length(args) >= 2 ? args[1:2] : ["35_01", "35_17"]
    stride = length(args) >= 3 ? parse(Int, args[3]) : 1
    trials = Dict(t.id => t for t in load_trials())
    a, b = map(ids) do id
        haskey(trials, id) || error("trial $id is not in trials.txt")
        trials[id]
    end
    return compare(a, b, Val(curve_dim(a)), stride)
end

main(ARGS) || exit(1)

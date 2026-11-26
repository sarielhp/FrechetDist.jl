#!/usr/bin/env julia
# All-pairs Fréchet distances between the motion-capture trials of trials.txt (walks,
# runs and jumps of five subjects), as curves of joint angles in ℝ^46.
#
# Prints, for each trial, its nearest other trial and its nearest trial of another
# subject, how often their activity agrees with the trial's own, and the mean distance
# between activities. Writes the full matrix to example_higher/output/distances.csv.
#
# Usage (from the repository root, after `ruby example_higher/download.rb`):
#   julia example_higher/distance_matrix.jl [STRIDE [APPROX]]
# STRIDE keeps every STRIDE-th frame of the 120 Hz recording (default 1). With APPROX
# (> 1, e.g. 1.05) the distances come from frechet_c_approx, within that factor.

using QuickEnv # local
using FrechetDist
using Printf

include(joinpath(@__DIR__, "mocap.jl"))
using .Mocap

const OUTPUT = joinpath(@__DIR__, "output", "distances.csv")

subject(t::Trial) = first(split(t.id, '_'))

function distance_matrix(trials, ::Val{D}, stride::Int, approx) where {D}
    curves = [trial_curve(t, Val(D); stride) for t in trials]
    n = length(curves)
    dist = zeros(n, n)
    for i in 1:n, j in (i+1):n
        m = if approx === nothing
            frechet_c_compute(curves[i], curves[j])
        else
            frechet_c_approx(curves[i], curves[j], approx)
        end
        dist[i, j] = dist[j, i] = m.leash
    end
    return dist
end

# The index of the trial nearest to trial i among those `admissible` says to consider.
function nearest(dist, i, admissible)
    candidates = [j for j in axes(dist, 2) if j != i && admissible(j)]
    return candidates[argmin(dist[i, candidates])]
end

function report_neighbours(trials, dist)
    println("trial  activity | nearest trial          | nearest of another subject")
    hits = [0, 0]
    for (i, t) in enumerate(trials)
        nn = nearest(dist, i, _ -> true)
        other = nearest(dist, i, j -> subject(trials[j]) != subject(t))
        hits .+= (trials[nn].activity == t.activity, trials[other].activity == t.activity)
        @printf(
            "%s  %-8s | %s %-5s %9.2f | %s %-5s %9.2f\n",
            t.id,
            t.activity,
            trials[nn].id,
            trials[nn].activity,
            dist[i, nn],
            trials[other].id,
            trials[other].activity,
            dist[i, other]
        )
    end
    n = length(trials)
    @printf(
        "\nSame activity as the nearest trial: %d/%d; as the nearest of another subject: %d/%d\n",
        hits[1],
        n,
        hits[2],
        n
    )
end

function report_activities(trials, dist)
    acts = unique(t.activity for t in trials)
    println("\nMean distance between activities:")
    @printf("%8s", "")
    foreach(a -> @printf("%10s", a), acts)
    println()
    for a in acts
        @printf("%-8s", a)
        for b in acts
            pairs = [
                dist[i, j] for i in eachindex(trials), j in eachindex(trials) if
                i != j && trials[i].activity == a && trials[j].activity == b
            ]
            @printf("%10.2f", sum(pairs) / length(pairs))
        end
        println()
    end
end

function write_csv(trials, dist)
    mkpath(dirname(OUTPUT))
    open(OUTPUT, "w") do io
        println(io, join(["trial"; [t.id for t in trials]], ','))
        for (i, t) in enumerate(trials)
            println(io, join([t.id; [@sprintf("%.6f", d) for d in dist[i, :]]], ','))
        end
    end
    println("\nDistance matrix written to $OUTPUT")
end

function main(args)
    stride = isempty(args) ? 1 : parse(Int, args[1])
    approx = length(args) >= 2 ? parse(Float64, args[2]) : nothing
    trials = load_trials()
    D = curve_dim(trials[1])
    method = approx === nothing ? "exact" : "$approx-approximate"
    println(
        "$(length(trials)) trials, curves in ℝ^$D, stride $stride, $method Fréchet distance",
    )
    seconds = @elapsed dist = distance_matrix(trials, Val(D), stride, approx)
    @printf("%d pairs in %.1f s\n\n", length(trials) * (length(trials) - 1) ÷ 2, seconds)
    report_neighbours(trials, dist)
    report_activities(trials, dist)
    write_csv(trials, dist)
end

main(ARGS)

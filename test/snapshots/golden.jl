# Usage: julia --project=<proj> golden.jl <out.jls>
# Dumps plain-number snapshots of the curve structures FrechetDist builds, so two
# revisions can be compared bit for bit (see compare.rb).
using FrechetDist, Random, Serialization
using BasicCompGeometry
import BasicCompGeometry: cardin, Points
const FD = FrechetDist
mk(D, pts) = PntSeq{D,Float64}([Point{D,Float64}(p...) for p in pts])
mat(P) = [Points(P)[i][d] for d = 1:length(Points(P)[1]), i = 1:cardin(P)]
conv(x::Union{FD.Morphing}) = (
    x.leash,
    mat(x.P),
    mat(x.Q),
    [(Vector(e.p), e.i, Int(e.type), e.t, e.penalty) for e in x.pes],
    [(Vector(e.p), e.i, Int(e.type), e.t, e.penalty) for e in x.qes],
)
conv(x::Tuple) = map(conv, x)
conv(x::AbstractVector{<:Real}) = collect(x)
conv(x::Real) = x
conv(x) = (hasproperty(x, :pnts) || applicable(cardin, x)) ? mat(x) : x
conv(x::FD.PolygonHierarchy) =
    (x.len, mat(x.P), collect(x.plt), [mat(q) for q in x.polys], collect(x.widths))

function walk(D, n, seed)
    rng = MersenneTwister(seed)
    cur = zeros(D)
    pts = Vector{Vector{Float64}}()
    for _ = 1:n
        cur = cur .+ (rand(rng, D) .- 0.5) .* (1 + 3 * rand(rng))
        push!(pts, copy(cur))
    end
    return mk(D, pts)
end

out = Dict{String,Any}()
rec(k, f) = (out[k] = conv(f()))
for (D, n, seed) in [(2, 3, 1), (2, 5, 2), (2, 60, 3), (2, 400, 4), (3, 120, 5), (4, 90, 6)]
    P = walk(D, n, seed)
    Q = walk(D, max(n - 7, 3), seed + 100)
    id = "D$D-n$n"
    rec("$id/width", () -> FD.frechet_width_approx(P))
    n > 10 && rec("$id/width_rng", () -> FD.frechet_width_approx(P, 3:(n÷2)))
    rec("$id/plt", () -> FD.frechet_palette(P))
    plt = FD.frechet_palette(P)
    for w in (0.0, 0.3, 1.0, 4.0, 1e9)
        rec("$id/plt_approx/$w", () -> FD.frechet_palette_approx(P, plt, w))
    end
    for lv in (0, 1, 3)
        rec("$id/plt_level/$lv", () -> redirect_stdout(devnull) do
            FD.frechet_palette_level(P, plt, lv)
        end)
    end
    for w in (0.2, 1.0, 5.0)
        rec("$id/simp_w_exp/$w", () -> FD.frechet_simplify_w_exp(P, w))
        rec("$id/simp_to_width/$w", () -> FD.frechet_simplify_to_width(P, w))
        rec("$id/simp_ext/$w", () -> FD.Polygon_simplify_ext(P, w))
        R, ri = FD.frechet_simplify_w_exp(P, w)
        rec("$id/mono_subcurve/$w", () -> FD.frechet_c_mono_approx_subcurve(P, R, ri))
        rec("$id/offsets/$w", () -> FD.frechet_offsets(P, ri))
    end
    rec(
        "$id/radii",
        () -> FD.Polygon_simplify_radii(P, [0.1 + 0.05 * (i % 7) for i = 1:cardin(P)]),
    )
    rec("$id/split_edges", () -> FD.split_edges(P))
    rec("$id/sample_unif", () -> FD.Polygon_sample_uniformly(P, 17))
    rec("$id/spine", () -> FD.spine(P))
    rec("$id/reverse", () -> Base.reverse(P))
    n >= 60 && rec("$id/simp_to_cardin", () -> redirect_stdout(devnull) do
        FD.frechet_simplify_to_cardin(P, 10)
    end)
    # Hierarchy
    if n >= 60
        ph = FD.ph_init(P, ReentrantLock())
        rec("$id/ph_init", () -> ph)
        FD.ph_compute_hierarchy(ph)
        rec("$id/ph_hier", () -> ph)
        for w in (0.5, 2.0, 8.0)
            rec("$id/ph_approx/$w", () -> FD.ph_approx(ph, w))
            rec("$id/ph_after/$w", () -> ph)
        end
    end
    # Distance algorithms (outputs of the main computations)
    rec("$id/c_compute", () -> FD.frechet_c_compute(P, Q))
    rec("$id/d_compute", () -> FD.frechet_d_compute(P, Q))
    rec("$id/ve_r", () -> FD.frechet_ve_r_compute(P, Q))
    rec("$id/sweep", () -> FD.SweepDist_compute(P, Q))
    rec("$id/dtw", () -> FD.DTW_d_compute(P, Q))
    mc = FD.frechet_c_compute(P, Q)
    rec("$id/sd_price", () -> FD.Morphing_SweepDist_price(mc))
    rec("$id/sd_approx_price", () -> FD.Morphing_SweepDist_approx_price(mc))
    rec("$id/sd_price_sweep", () -> FD.Morphing_SweepDist_price(FD.SweepDist_compute(P, Q)))
end

# File readers
mktempdir() do dir
    files = Dict(
        "comma" => "1.5,2.5\n1.5,2.5\n3,4\n5.25,-6\n",
        "space" => "1 2\n3 4\n3 4\n7 9\n",
        "header" => "x y k tid\n1 2 3 4\n5 6 7 8\n",
        "single" => "4,5\n",
        "dups" => "1,1\n1,1\n",
    )
    for (nm, body) in files
        f = joinpath(dir, nm)
        write(f, body)
        rec("read_file/$nm", () -> FD.read_file(f))
    end
    g = joinpath(dir, "geo.plt")
    write(
        g,
        "h1\nh2\nh3\nh4\nh5\nh6\n39.9,116.3,0,492,39744.12,2008-10-23,02:53:04\n39.9,116.3,0,492,39744.12,2008-10-23,02:53:05\n39.95,116.4,0,492,39744.13,2008-10-23,02:53:06\n",
    )
    rec("read_geolife", () -> FD.read_geolife_plt(g))
end
serialize(ARGS[1], out)
println("wrote ", length(out), " entries")

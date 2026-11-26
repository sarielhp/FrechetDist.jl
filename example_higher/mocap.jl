# Reader of the CMU motion-capture trials fetched by download.rb. A trial (Acclaim .amc
# file) lists, frame by frame, the joint angles of a skeleton in degrees; each frame
# becomes one point of a curve whose dimension is the number of angles.
#
# Some bones are dropped. The root: its six values are the position in the lab and the
# global heading, which say where the subject was and which way it faced, not how it
# moved. The hands, fingers and thumbs: the database tracks them with few markers, and
# their angles show Euler flips and single-frame spikes of over 150 degrees. Every
# other bone keeps its angles, ordered by bone name, so curves of different trials and
# subjects share coordinates whenever their skeletons have the same bones.
#
# Missing data is written as all-zero values of a bone (16_15 starts with two such
# frames for the left humerus). The Fréchet distance is a bottleneck measure, so a
# single such frame can decide it; frames where a bone that is usually nonzero reads
# exactly zero are dropped.

module Mocap

using BasicCompGeometry: Point, PntSeq

export Trial, load_trials, trial_curve, curve_dim

const DATA_DIR = joinpath(@__DIR__, "data")
const MANIFEST = joinpath(@__DIR__, "trials.txt")
const DROPPED_BONES = ("root", "lhand", "rhand", "lfingers", "rfingers", "lthumb", "rthumb")

"""
    Trial

One motion: its id (`"35_01"`), activity class from `trials.txt`, the bones with their
number of angles, the angles as a matrix with one column per frame, and the number of
frames dropped for missing data.
"""
struct Trial
    id::String
    activity::String
    layout::Vector{Pair{String,Int}}
    angles::Matrix{Float64}
    dropped::Int
end

Base.length(t::Trial) = size(t.angles, 2)
curve_dim(t::Trial) = size(t.angles, 1)

"""
    read_amc_frames(path) -> Vector{Dict{String,Vector{Float64}}}

The frames of an .amc file, each a map from bone name to its values.
"""
function read_amc_frames(path::AbstractString)
    frames = Dict{String,Vector{Float64}}[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || line[1] in ('#', ':')) && continue
        if all(isdigit, line)
            push!(frames, Dict{String,Vector{Float64}}())
            continue
        end
        isempty(frames) && error("$path: bone values before the first frame number")
        fields = split(line)
        frames[end][fields[1]] = parse.(Float64, fields[2:end])
    end
    isempty(frames) && error("$path: no frames")
    return frames
end

function frame_layout(frame::Dict{String,Vector{Float64}})
    bones = sort!([b for b in keys(frame) if !(b in DROPPED_BONES)])
    return [b => length(frame[b]) for b in bones]
end

"""
    read_trial(id, activity) -> Trial

Reads `data/<id>.amc`. Every frame must carry the same bones with the same number of
angles as the first one.
"""
function read_trial(id::AbstractString, activity::AbstractString)
    path = joinpath(DATA_DIR, "$id.amc")
    isfile(path) || error("$path is missing; run `ruby example_higher/download.rb` first")
    frames = read_amc_frames(path)
    layout = frame_layout(frames[1])
    angles = Matrix{Float64}(undef, sum(last, layout), length(frames))
    for (j, frame) in enumerate(frames)
        frame_layout(frame) == layout || error("$path: frame $j has another set of bones")
        row = 0
        for (bone, k) in layout
            angles[(row+1):(row+k), j] = frame[bone]
            row += k
        end
    end
    keep = complete_frames(layout, angles)
    any(keep) || error("$path: every frame has missing data")
    return Trial(id, activity, layout, angles[:, keep], count(!, keep))
end

"""
    complete_frames(layout, angles) -> BitVector

The frames without missing data: a frame is missing a bone when all the bone's angles
are exactly zero there although they are not in most frames of the trial.
"""
function complete_frames(layout::Vector{Pair{String,Int}}, angles::Matrix{Float64})
    keep = trues(size(angles, 2))
    row = 0
    for (_, k) in layout
        rows = (row+1):(row+k)
        row += k
        zero_frames = [all(iszero, view(angles, rows, j)) for j in axes(angles, 2)]
        2 * count(zero_frames) < length(zero_frames) || continue
        keep .&= .!zero_frames
    end
    return keep
end

"""
    load_trials(; activities = nothing) -> Vector{Trial}

The trials of `trials.txt`, optionally only those whose activity is in `activities`.
All of them must have the same bone layout, hence curves of the same dimension.
"""
function load_trials(; activities = nothing)
    trials = Trial[]
    for raw in eachline(MANIFEST)
        fields = split(replace(raw, r"#.*" => ""))
        isempty(fields) && continue
        id, activity = fields
        (activities === nothing || activity in activities) || continue
        push!(trials, read_trial(id, activity))
    end
    isempty(trials) && error("no trials selected from $MANIFEST")
    for t in trials
        t.layout == trials[1].layout ||
            error("trial $(t.id) has another skeleton than $(trials[1].id)")
    end
    return trials
end

"""
    trial_curve(t, Val(D); stride = 1) -> PntSeq{D,Float64}

The trial as a polygonal curve in ``ℝ^D``, keeping every `stride`-th frame (and the
last one). `D` must equal `curve_dim(t)`; passing it as a `Val` makes the curve type,
and so the distance computations, type-stable.
"""
function trial_curve(t::Trial, ::Val{D}; stride::Int = 1) where {D}
    D == curve_dim(t) || error("trial $(t.id) has dimension $(curve_dim(t)), not $D")
    cols = collect(1:stride:length(t))
    last(cols) == length(t) || push!(cols, length(t))
    P = PntSeq{D,Float64}()
    for j in cols
        push!(P, Point{D,Float64}(ntuple(i -> t.angles[i, j], Val(D))))
    end
    return P
end

end # module Mocap

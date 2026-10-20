using FrechetDist
using BasicCompGeometry
using Test

# Run `code` in a fresh Julia process on this project; it must exit with status 0 within
# `seconds`. A computation that never terminates is killed, so it fails the test instead of
# hanging the run (or starving the test's own threads).
function finishes(code::AbstractString, seconds)
    project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --project=$project -e $code`
    proc = run(pipeline(cmd; stdout = devnull, stderr = devnull); wait = false)
    deadline = time() + seconds
    while process_running(proc) && time() < deadline
        sleep(0.2)
    end
    if process_running(proc)
        kill(proc, Base.SIGKILL)
        wait(proc)
        return false
    end
    return success(proc)
end

const PRELUDE = "using FrechetDist, BasicCompGeometry; const FD = FrechetDist; "

@testset "frechet_simplify_to_cardin terminates on zero-width curves" begin
    collinear = PRELUDE * """
        P = PntSeq([point(Float64(i), 0.0) for i = 1:60])
        Q, idx = FD.frechet_simplify_to_cardin(P, 10)
        exit(length(Q.pnts) == 2 && idx == [1, 60] ? 0 : 3)
        """
    @test finishes(collinear, 60)

    same = PRELUDE * """
        FD.frechet_simplify_to_cardin(PntSeq([point(1.0, 1.0) for _ = 1:20]), 3)
        """
    @test finishes(same, 60)

    # Three corners joined by long collinear runs: the result keeps every corner and is exact
    # (it may keep a few redundant vertices on the straight runs).
    corners = PRELUDE * """
        runs = [point(Float64(i), 0.0) for i = 0:30]
        append!(runs, [point(30.0, Float64(i)) for i = 1:30])
        append!(runs, [point(30.0 - i, 30.0) for i = 1:30])
        Q, idx = FD.frechet_simplify_to_cardin(PntSeq(runs), 10)
        corners = [1, 31, 61, 91]
        exit(all(c -> c in idx, corners) && length(Q.pnts) < 12 ? 0 : 3)
        """
    @test finishes(corners, 60)

    # A curve that does need the refinement loop still returns a simplification.
    spiral = PRELUDE * """
        P = PntSeq([point(cos(0.3i) * i, sin(0.3i) * i) for i = 1:60])
        Q, idx = FD.frechet_simplify_to_cardin(P, 10)
        exit(length(Q.pnts) >= 12 || 2 * length(Q.pnts) >= 60 ? 0 : 3)
        """
    @test finishes(spiral, 120)
end

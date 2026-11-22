# Temporary files of the examples live in <tmpdir>/<user>/FrechetDist/<purpose>/,
# per the standards' temporary-file policy. This is the only place that calls
# tempdir().
function example_tmp(purpose::AbstractString)
    user = get(ENV, "USER", string(Base.Libc.getuid()))
    base = joinpath(tempdir(), user)  # tmp-exempt: the helper
    mkpath(base; mode = 0o700)
    dir = joinpath(base, "FrechetDist", purpose)
    mkpath(dir; mode = 0o700)
    return dir
end

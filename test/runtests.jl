using FrechetDist
using BasicCompGeometry
using Test

function example_01()
    P = PntSeq2F()
    Q = PntSeq2F()

    push!(P, point(0.0, 0), point(1.0, 0.0))
    push!(Q, point(0.0, 1.0), point(1.0, 1.0))
    return P, Q
end


function example_02()
    P = rand_pnt_seq(3, Float64, 100)
    Q = rand_pnt_seq(3, Float64, 80)

    return P, Q
end

#@testset "FrechetDist.jl" begin
@testset "Frechet Dist.jl" begin

    # Write your tests here.
    P, Q = example_01()

    @test FrechetDist.frechet_c_compute(P, Q).leash == 1.0

    P, Q = example_02()

    # Code generation - forcing to generate all different distance
    # calculations
    m_d = frechet_d_compute(P, Q)
    @test m_d.leash > 0.0

    m_d_r = frechet_d_r_compute(P, Q)
    @test m_d_r.leash > 0.0
    @test m_d.leash == m_d_r.leash

    m_ve = frechet_ve_r_compute(P, Q)
    @test m_ve.leash > 0.0

    m_e = frechet_c_compute(P, Q)
    @test m_e.leash > 0.0

    @test (1.000001 * m_e.leash) >= m_ve.leash

    m_a_2 = frechet_c_approx(P, Q, 2.0)

    @test (m_a_2.leash > 0.0)

    @test ((m_e.leash <= m_a_2.leash) && (m_a_2.leash <= 2.0*m_e.leash))

    d = frechet_dist_upper_bound(P, Q)
    @test m_e.leash <= d

    m_mono, f, P_r, Q_r = frechet_mono_via_refinement(P, Q, 1.01)
    @test 1.000001*m_mono.leash >= m_ve.leash

    w_P = frechet_width_approx(P)

    w_Q = frechet_width_approx(Q)
    @test w_P >= 0
    @test w_Q >= 0
end

# Regression: frechet_c_compute used to return a value up to 1e-3 too large on long,
# nearly collinear GPS tracks (it stopped refining too early). The curves are the first
# 200 GPS points of two GeoLife trajectories, in meters. A prefix with n and m vertices
# corresponds to the first 120/160/180/200 raw GPS points of each file.
@testset "frechet_c_compute on GPS tracks" begin
    P = FrechetDist.read_file(joinpath(@__DIR__, "fixtures", "gps_P.txt"))
    Q = FrechetDist.read_file(joinpath(@__DIR__, "fixtures", "gps_Q.txt"))
    @test (length(P), length(Q)) == (102, 200)

    for (n, m) in ((57, 120), (68, 160), (88, 180), (102, 200))
        Pn = PntSeq(P.pnts[1:n])
        Qm = PntSeq(Q.pnts[1:m])
        # An upper bound within 1e-7 (relative) of the true distance.
        a = frechet_c_approx(Pn, Qm, 1.0000001).leash
        for f_accept_approx in (false, true)
            c = frechet_c_compute(Pn, Qm, f_accept_approx).leash
            @test c >= a * (1 - 1e-6)
            @test c <= a * (1 + 2e-5)
        end
    end
end

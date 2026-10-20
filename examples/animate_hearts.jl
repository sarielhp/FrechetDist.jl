#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist, Plots
using BasicCompGeometry

# Appends the points f( r ), for r in range, to P.
function fill_curve(P::PntSeq{D,T}, f, range) where {D,T}
    for r in range
        push!(P, Point{D,T}(f(r)...))
    end
    return P
end

# t ∈ [-1,1]
#heart(t) = 

#P = fill_curve( Polygon2F(), x->( sin(x), cos(x) )    , 0:0.01:2*pi )
#Q = fill_curve( Polygon2F(), x->( sin(2.0*x), cos(x) ), 0:0.01:2*pi )


# t ∈ [-1, 1]
heart(t) = (3.0 * sin(t) * cos(t) * log(abs(t)), 3.0*abs(t)^0.3*(cos(t))^0.5 - 1.5)


uheart(t) = (3.0 * sin(t) * cos(t) * log(abs(t)), -3.0*abs(t)^0.3*(cos(t))^0.5 + 1.5)

heart(t) = (-16*(sin(t))^3, 13.0*cos(t)-5.0*cos(2*t) - 2.0 * cos(3*t) - cos(4.0*t))
uheart(t) = (16*(sin(t))^3, -(13.0*cos(t)-5.0*cos(2*t) - 2.0 * cos(3*t) - cos(4.0*t)))

P = fill_curve(PntSeq2F(), heart, (-pi):0.02:pi);
Q = fill_curve(PntSeq2F(), uheart, 0:0.02:(2.0*pi));

mrp = frechet_c_compute(P, Q);

M = Morphing_sample_uniformly(mrp, 200);

anim = @animate for i = 1:size(M, 1)
    default(
        grid = :false,
        legend = false,
        framestyle = :none,
        label = :none,
        ticks = false,
        showaxis = false,
        lw = 4,
    )
    plt = plot(300, 300, 0, dpi = 200)
    plot!([p[1] for p in P], [p[2] for p in P], lc = :red)
    plot!([p[1] for p in Q], [p[2] for p in Q], lc = :green)

    t = [M[i, 1:2]'; M[i, 3:4]']
    plot!(t[:, 1], t[:, 2], lw = 8, lc = :blue)
end

gif(anim, "frechet.gif", fps = 10)

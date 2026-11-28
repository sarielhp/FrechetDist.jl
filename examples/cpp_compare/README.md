# FrechetDist against the C++ code of Bringmann, Künnemann and Nusser

The decider benchmark of K. Bringmann, M. Künnemann and A. Nusser, *Walking the Dog
Fast in Practice: Algorithm Engineering of the Fréchet Distance* (SoCG 2019), run on
their C++ implementation
([gitlab.mpi-klsb.mpg.de/anusser/frechet_distance](https://gitlab.mpi-klsb.mpg.de/anusser/frechet_distance),
commit `93f8d44`) and on FrechetDist, on the same curve files and the same queries.

Each query is a pair of curves and a distance `r`; the answer is whether the Fréchet
distance of the curves is below `r`. The queries come in files
`<set>_query_decider_<k>_<l>_<plus|minus>.txt` of 1000 queries each, whose distance is
the Fréchet distance of the pair times `1 + 10^l` (plus) or `1 - 10^l` (minus): the
smaller `l`, the closer the query is to the distance and the harder it is to decide.
There are 897 files: 322 for sigspatial (20,199 curves), 253 for characters (2,858
curves) and 322 for geolife (16,966 curves).

## Results: exact decisions

One thread, AMD Ryzen 7 5700G, Julia 1.13.3, g++ 16.2. Both programs read the curves
of a query file first, with the per-curve work their deciders rely on (the C++ `Curve`
computes the bounding box and prefix lengths of a curve; FrechetDist's
`DecideCurve(P, BoxTree)` checks the curve and computes its width, bounding box, prefix
lengths and a tree of the bounding boxes of its pieces), run the queries once
untimed (warm-up), and then time the decision calls only: the minimum of three runs
per query file. Every program gives the same answer as the C++ code to all 897,000
queries. The FrechetDist column was measured with trees down to single segments and
4 x 4 leaf boxes (the trees now stop at runs of 8 segments, which is a little faster;
see below).

| Dataset | Queries | C++ | `frechet_c_decide_exact` | ratio | `run_queries_test.jl` | ratio |
|:--|--:|--:|--:|--:|--:|--:|
| sigspatial | 322,000 | 3.4 s | 2.1 s | 0.6 | 100.1 s | 29 |
| characters | 253,000 | 5.2 s | 1.8 s | 0.3 | 72.2 s | 14 |
| geolife | 322,000 | 18.9 s | 16.5 s | 0.9 | 1414.6 s | 75 |

Reading the curves (not timed above) takes 36.8 s, 21.4 s and 155.7 s for C++, and
51.0 s, 23.0 s and 200.0 s for FrechetDist, which also builds the box trees of the
curves (each query file reads its curves again).

`frechet_c_decide_exact` is the library's exact decision procedure, on prepared curves
(`FD_MODE=exact`; `exact_lazy` prepares them with the default `DecideCurve(P)`, see
below, and `exact_plain` calls it on the curves themselves);
`run_queries_test.jl` is the decider `frechet_decider_PID` of
`examples/run_queries_test.jl` (`FD_MODE=pid`, the default; measured before the
warm-up pass was added, the minimum of three runs, with its simplification
hierarchies built when the curves are read).

`frechet_c_decide_exact` tries, in order: the distances between the endpoints; the
largest distance between the bounding boxes of the curves; the endpoint distance plus
the widths of the curves; a greedy monotone walk through pairs of vertices, and one by
adaptive steps over pieces of the curves that lie in balls within `r` of each other
(both prove "yes"); a vertex of one curve with no point of the other within `r`
(proves "no"); decisions on simplifications of the curves (error `w`, from `r/8`
down, while a simplification keeps at most an eighth of the vertices), which settle
the query when the simplified curves are more than `2w` from `r`; and the free space
of the curves themselves. The free space is searched by a recursive decomposition
into boxes (in four when square, otherwise in two across the longer side): a box whose
two curve pieces are all within `r` of each other (or all
farther apart) is settled through its boundary, a box that nothing enters is skipped,
and boxes of at most 8 x 8 cells are visited cell by cell. All the tests are exact, up
to the rounding of the free intervals.

Per group of query files, FrechetDist takes 0.3 to 1.0 times the time of the C++
code. The groups where the two are closest are the geolife queries just above the
distance (`l` from -5 to -10), where nearly all of the time goes to the free-space
search.
The breakdown by `l` is in
`raw/cpp_compare/summary_dexact.txt` after a run
(`PROGRAMS=fd FD_MODE=exact TAG=_dexact ruby examples/cpp_compare/run.rb 3`, then
`ruby examples/cpp_compare/summarize.rb _dexact`).

### Trees built inside the timed decisions

The default `DecideCurve(P)` builds the tree of a curve on demand (`LazyBoxTree`): the
first read builds its top (nodes of at least 128 segments) by one scan of the curve,
and the first read of a smaller node builds its subtree, down to runs of 8 segments, by
one scan of its piece; the prefix lengths are computed when first needed. So the tree
building is part of the timed decisions (`FD_MODE=exact_lazy`: the curves are prepared
again, untimed, before each timed run), and the C++ code still computes its bounding
boxes and prefix lengths untimed. On 40 random query files per dataset (all `k` and
`l`):

| Dataset | C++ | trees built when read (untimed) | + their build time | trees built on demand (timed) | ratio to C++ |
|:--|--:|--:|--:|--:|--:|
| sigspatial | 0.260 s | 0.163 s | 0.125 s | 0.176 s | 0.68 |
| characters | 0.388 s | 0.141 s | 0.042 s | 0.148 s | 0.38 |
| geolife | 1.575 s | 1.234 s | 0.396 s | 1.342 s | 0.85 |

Most decisions are settled by the filters or near the root of the trees, so building
on demand does 10% to 27% of the work of building every tree when the curves are read.
FrechetDist is slower than the C++ code on the easier geolife files (`k` = 14, up to
2.4 times on single files), where the free-space search itself costs more.

## Setup

From the root of the repository (`raw/` is not tracked):

1. Clone the C++ code into `raw/frechet` and fetch the data with their URLs
   (`test_data/benchmark/fetch_and_convert_data.sh`): the GeoLife zip, `shortest-sf.tgz`
   and `mixoutALL_shifted.mat`.
2. Convert with their scripts, unchanged: `geolife_converter.py` (in
   `Geolife Trajectories 1.3/`), the sigspatial loop of `fetch_and_convert_data.sh`, and
   `character_converter.m`. Without octave, the characters files can be made with a
   port of the octave script that writes the same bytes (each row is the sum of the
   previous velocity rows, written as `" %.8e %.8e"`).
3. Build the driver with the flags of their CMake default build (RelWithDebInfo), from
   `raw/frechet`:

   ```
   g++ -O2 -g -DNDEBUG -std=c++11 -fno-omit-frame-pointer -DNVERBOSE -fopenmp \
       -DWITH_OPENMP -flto -Isrc -o build_bench/decider_bench \
       ../../examples/cpp_compare/decider_bench.cpp src/frechet_light.cpp \
       src/frechet_naive.cpp src/geometry_basics.cpp src/filter.cpp \
       src/orth_range_search.cpp src/parser.cpp src/query.cpp src/times.cpp src/curve.cpp
   ```

4. `ruby examples/cpp_compare/run.rb 3`, then `ruby examples/cpp_compare/summarize.rb`
   (environment variables `FD_MODE`, `PROGRAMS` and `TAG`: see `run.rb`).

`decider_bench.cpp` times `FrechetLight::lessThanWithFilters` as their
`deciderComparisonExp` does (curves read first, one `FrechetLight` object for all
queries), except that each query file is run once untimed before the timed runs, and
the minimum of the runs is reported (theirs averages five runs). `fd_bench.jl` runs
`frechet_c_decide_exact` on `DecideCurve`s, or the decider of
`examples/run_queries_test.jl` (`frechet_decider_PID`, with the hierarchies of its
serial mode), on the same files, in the same way.
Both programs read the first two numbers of each line and drop repeated consecutive
points, so they see the same curves.

## Pitfalls

- **Geolife numbering.** The geolife queries name curves by number (`5228.txt`), and
  the numbers come from the order in which `geolife_converter.py` writes the
  trajectories near Beijing; trajectories that leave the area are dropped without a
  number. The modified converter in `examples/replacement/` keeps those trajectories
  under other names, so every later number refers to a different trajectory (18,556
  of the files differ). Results on data converted that way are meaningless: use the
  original converter.
- **No `-march=native` for the C++ code.** Built with `-march=native`, g++ contracts
  multiplications and additions into FMA instructions, and the C++ deciders answer 17
  geolife queries wrongly (distance above `r` by 0.1% to 3%, confirmed by FrechetDist's
  exact algorithms and by the C++ code's own naive decider). With
  `-ffp-contract=off`, or with their default flags, the answers are right.

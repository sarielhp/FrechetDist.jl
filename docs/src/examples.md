# Examples

Here we describe how to download test date, and run some of the
examples. Here, I assume you are installing the package directly from
github, and not from Julia package manager.

## Getting started

To get the source code, do
```
mkdir FrechetDist
cd FrechetDist
git clone git@github.com:sarielhp/FrechetDist.jl.git .
```

Assuming the Julia is already installed, to install prerequisite packages,
do:
```
julia examples/prereq.jl
```
(This would take a few minutes potentially.)

## Getting the data

Do:

```julia examples/download_data.jl```

This creates `data/` subdirectory with a lot of data from various
sources (this directory is almost 1GB).  (This would take a few
minutes potentially.) This also creates a `raw/` subdirectory, which
you absolutely should delete:

```rm -r -f raw/```

Since it is (3.6GB), and contains all the downloaded data before conversion.

## Testing the data

We use test data available from earlier
[project](https://gitlab.mpi-klsb.mpg.de/anusser/frechet_distance),
which in turn was inspired (or is from) sigspatial competition. In any
case, you have three scripts to run the tests, the larger one takes
several hours:
```
examples/scripts/run_test_queries_characters
examples/scripts/run_test_queries_geolife
examples/scripts/run_test_queries_sigspatial
```
If you change the string `sfile` in these scripts to `file`, the
scripts would perform the tests in parallel using all available
threads, and would be ***much*** faster.

## Generating the examples in the webpage

To get more interesting outputs, and general the example webpage of
this page (see
[here](https://sarielhp.org/p/24/frechet_ve/examples/)), you can run
the command:
```
julia examples/generate_examples.jl
```

## For a minimal example...

See the `test/` subdirectory.

#!/usr/bin/env ruby
# frozen_string_literal: true

# Runs the decider benchmark of Bringmann, Künnemann and Nusser on their C++ code and on
# FrechetDist, one program at a time on one thread. See README.md for the setup.
#
#   ruby examples/cpp_compare/run.rb [runs] [dataset...]
#
# Outputs go to raw/cpp_compare/: <dataset>_<program><TAG>.log (one line per query
# file) and <dataset>_<program><TAG>_answers.txt (one line per query). Environment:
# FD_MODE is the mode of fd_bench.jl (default pid), TAG a suffix for the output names
# (default none), PROGRAMS the programs to run (default cpp,fd).

ROOT = File.expand_path('../..', __dir__)
CPP = File.join(ROOT, 'raw/frechet')
BENCH = File.join(CPP, 'test_data/benchmark')
QUERIES = File.join(CPP, 'test_data/decider_benchmark_queries')
OUT = File.join(ROOT, 'raw/cpp_compare')
CURVES = {
  'sigspatial' => File.join(BENCH, 'sigspatial/'),
  'characters' => File.join(BENCH, 'characters/data/'),
  'geolife' => File.join(BENCH, 'Geolife Trajectories 1.3/data/')
}.freeze

runs = ARGV.first || '3'
sets = ARGV.drop(1)
sets = CURVES.keys if sets.empty?
tag = ENV.fetch('TAG', '')
fd_mode = ENV.fetch('FD_MODE', 'pid')
programs = ENV.fetch('PROGRAMS', 'cpp,fd').split(',')
Dir.mkdir(OUT) unless Dir.exist?(OUT)
env = { 'JULIA_NUM_THREADS' => '1', 'OMP_NUM_THREADS' => '1' }

sets.each do |set|
  files = Dir.glob(File.join(QUERIES, "#{set}_query_decider_*.txt")).sort
  answers = ->(prog) { File.join(OUT, "#{set}_#{prog}#{tag}_answers.txt") }
  cmds = {
    'cpp' => [File.join(CPP, 'build_bench/decider_bench'), CURVES[set], answers['cpp'], runs, *files],
    'fd' => ['julia', "--project=#{ROOT}", File.join(__dir__, 'fd_bench.jl'), fd_mode, CURVES[set],
             answers['fd'], runs, *files]
  }
  cmds.slice(*programs).each do |prog, cmd|
    warn "#{set}: #{prog} #{prog == 'fd' ? fd_mode : ''} (#{files.size} query files)"
    log = File.join(OUT, "#{set}_#{prog}#{tag}.log")
    system(env, *cmd, out: log) or abort "#{set}: #{prog} failed"
  end
end

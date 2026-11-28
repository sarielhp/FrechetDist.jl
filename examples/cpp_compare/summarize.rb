#!/usr/bin/env ruby
# frozen_string_literal: true

# Summarizes the outputs of run.rb (raw/cpp_compare/): total decision time per dataset
# for each program, broken down by the query offset l of the file names
# (<set>_query_decider_<k>_<l>_<plus|minus>.txt: the query distance is the Fréchet
# distance times 1 + 10^l (plus) or 1 - 10^l (minus)), and every query on which the two
# programs disagree.
#
#   ruby examples/cpp_compare/summarize.rb [fd_tag eps]
#
# With a tag, the FrechetDist outputs of that run (TAG of run.rb) are compared with the
# untagged C++ outputs. With eps > 0 a disagreement is allowed in the files whose query
# distance is within a factor 1 +- eps of the distance (10^l <= eps), and counted
# separately; any other disagreement is an error.

OUT = File.expand_path('../../raw/cpp_compare', __dir__)
FD_TAG = ARGV[0] || ''
EPS = (ARGV[1] || '0').to_f

# file => [queries, yes, ms]; the loading times (ms, untimed) summed separately.
def read_log(path)
  rows = {}
  prep = 0.0
  File.foreach(path) do |line|
    f = line.split
    next unless f.size >= 4 && f[0].end_with?('.txt')

    rows[File.basename(f[0])] = [f[1].to_i, f[2].to_i, f[3].to_f]
    prep += f[4].to_f if f[4]
  end
  [rows, prep]
end

def read_answers(path)
  File.foreach(path).map do |line|
    a, b, r, ans = line.split
    [a, b, r.to_f, ans.to_i]
  end
end

def offset(file) = file[/_(-?\d+)_(plus|minus)\.txt\z/, 1].to_i
def sign(file) = file[/_(plus|minus)\.txt\z/, 1]

%w[sigspatial characters geolife].each do |set|
  cpp_log = File.join(OUT, "#{set}_cpp.log")
  fd_log = File.join(OUT, "#{set}_fd#{FD_TAG}.log")
  next unless File.exist?(cpp_log) && File.exist?(fd_log)

  cpp, cpp_prep = read_log(cpp_log)
  fd, prep = read_log(fd_log)
  common = cpp.keys & fd.keys
  n = common.sum { |f| cpp[f][0] }
  t_cpp = common.sum { |f| cpp[f][2] } / 1000
  t_fd = common.sum { |f| fd[f][2] } / 1000
  puts "== #{set}: #{common.size} query files, #{n} queries"
  puts format('  C++ %.2f s (+ %.2f s loading)   FrechetDist %.2f s (+ %.2f s loading)   ' \
              'FrechetDist / C++ = %.1f', t_cpp, cpp_prep / 1000, t_fd, prep / 1000, t_fd / t_cpp)
  common.group_by { |f| [sign(f), offset(f)] }.sort_by { |(s, l), _| [s == 'minus' ? 0 : 1, s == 'minus' ? -l : l] }
        .each do |(s, l), fs|
    c = fs.sum { |f| cpp[f][2] } / 1000
    d = fs.sum { |f| fd[f][2] } / 1000
    puts format('    %-5s l=%3d  C++ %8.3f s  FD %8.3f s  ratio %6.1f', s, l, c, d, d / c)
  end

  a = read_answers(File.join(OUT, "#{set}_cpp_answers.txt"))
  b = read_answers(File.join(OUT, "#{set}_fd#{FD_TAG}_answers.txt"))
  abort "#{set}: answer files have different queries" unless a.size == b.size &&
                                                          a.zip(b).all? { |x, y| x[0..2] == y[0..2] }
  # The query file of answer k: the logs list the files in order, with their sizes.
  files = File.foreach(cpp_log).flat_map { |l| f = l.split; [File.basename(f[0])] * f[1].to_i }
  bad = a.each_index.reject { |k| a[k][3] == b[k][3] }
  allowed, wrong = bad.partition { |k| EPS.positive? && 10.0**offset(files[k]) <= EPS * (1 + 1e-9) }
  puts "  disagreements: #{bad.size} of #{a.size}" +
       (EPS.positive? ? " (#{allowed.size} within the tolerance, #{wrong.size} outside it)" : '')
  allowed.map { |k| sign(files[k]) + ' l=' + offset(files[k]).to_s }.tally.sort.each do |key, n|
    puts "    within tolerance: #{key}: #{n}"
  end
  wrong.first(20).each do |k|
    puts "    #{files[k]}: #{a[k][0]} #{a[k][1]} r=#{a[k][2]}  C++ #{a[k][3]}  FD #{b[k][3]}"
  end
end

#!/usr/bin/env ruby
# frozen_string_literal: true

# Downloads the UCI PEMS-SF data (https://archive.ics.uci.edu/dataset/204/pems+sf, CC BY
# 4.0; 109 MB) into data/pems/, unless present, and converts it into curves: day k is a
# curve with 144 vertices (one per 10 minutes), each a point in R^963 (the occupancy
# rates, in [0, 1], of 963 San Francisco Bay Area freeway sensors). The 440 days span
# 2008-01-01 to 2009-03-30, without public holidays and two days of sensor failure;
# days are numbered in calendar order (randperm undoes the shuffle of the train/test
# files).
#
# Writes data/curves/pems/day<k>.f64 (963 x 144 Float64, column-major, little-endian:
# vertex after vertex) and data/curves/pems/index.csv (name,label,dims,vertices,file),
# the format read by curves.jl and fred_curves.py. The label is the weekday (1 =
# Monday ... 7 = Sunday).
#
# Usage: ruby example_higher/prepare_pems.rb

require 'fileutils'
require 'net/http'
require 'uri'

URL = 'https://archive.ics.uci.edu/static/public/204/pems+sf.zip'
SRC = File.join(__dir__, 'data', 'pems')
OUT = File.join(__dir__, 'data', 'curves', 'pems')
SENSORS = 963
STEPS = 144
DAY_NAMES = %w[Mon Tue Wed Thu Fri Sat Sun].freeze

# Streams the body of `url` into `target` through a partial file.
def download(url, target)
  partial = "#{target}.part"
  Net::HTTP.get_response(URI(url)) do |response|
    raise "#{url}: HTTP #{response.code}" unless response.is_a?(Net::HTTPSuccess)

    File.open(partial, 'wb') { |io| response.read_body { |chunk| io.write(chunk) } }
  end
  File.rename(partial, target)
end

def fetch_source
  return if File.exist?(File.join(SRC, 'randperm'))

  FileUtils.mkdir_p(SRC)
  zip = File.join(SRC, 'pems_sf.zip')
  unless File.exist?(zip)
    puts "downloading #{URL}"
    download(URL, zip)
  end
  system('unzip', '-o', '-q', zip, '-d', SRC, exception: true)
end

def read_vector(name)
  File.read(File.join(SRC, name)).delete('[]').split.map(&:to_i)
end

# One line: "[s1_t1 ... s1_t144;s2_t1 ...;...]", sensor rows. Returns the values in
# vertex order (all sensors at t1, then t2, ...), as Float64 bytes.
def day_bytes(line)
  rows = line.strip.delete('[]').split(';').map { |row| row.split.map(&:to_f) }
  unless rows.size == SENSORS && rows.all? { |r| r.size == STEPS }
    abort "#{SRC}: expected #{SENSORS} x #{STEPS} values per day"
  end
  rows.transpose.flatten.pack('E*')
end

def each_day(&)
  return enum_for(:each_day) unless block_given?

  %w[PEMS_train PEMS_test].each do |name|
    File.foreach(File.join(SRC, name), &)
  end
end

fetch_source
order = read_vector('randperm') # line k of train+test is calendar day order[k]
labels = read_vector('PEMS_trainlabels') + read_vector('PEMS_testlabels')
abort 'randperm and labels disagree in length' unless order.size == labels.size

FileUtils.mkdir_p(OUT)
rows = []
each_day.with_index do |line, k|
  day = order[k]
  file = format('day%03d.f64', day)
  File.binwrite(File.join(OUT, file), day_bytes(line))
  rows << [day, "day#{day}_#{DAY_NAMES[labels[k] - 1]}", labels[k], file]
end
abort "expected #{order.size} days, read #{rows.size}" unless rows.size == order.size

File.open(File.join(OUT, 'index.csv'), 'w') do |io|
  io.puts 'name,label,dims,vertices,file'
  rows.sort.each { |_, name, label, file| io.puts [name, label, SENSORS, STEPS, file].join(',') }
end
puts "#{rows.size} days written to #{OUT}"

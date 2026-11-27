#!/usr/bin/env ruby
# frozen_string_literal: true

# Downloads part of the KTH human action video database (Schuldt, Laptev and Caputo,
# ICPR 2004; https://www.csc.kth.se/cvap/actions/) into data/kth/, unless present, and
# turns videos into curves: every frame (160 x 120 pixels, grayscale, 25 fps) is a point
# in R^19200 (its pixels, scaled to [0, 1], row by row), and a curve is one repetition
# of an action, its frames in order (the first repetition listed in 00sequences.txt).
# Takes condition d1 (outdoors, static background) of persons 1 to 25 for the actions
# given (default walking and boxing; each archive is 150-240 MB). Needs ffmpeg.
#
# Writes data/curves/kth/<person>_<action>.f64 (19200 x frames Float64, column-major,
# little-endian: frame after frame) and data/curves/kth/index.csv
# (name,label,dims,vertices,file), the format read by curves.jl and fred_curves.py,
# persons first, so that curves 1..2k are k persons doing each action.
#
# Usage: ruby example_higher/prepare_kth.rb [ACTION ...]

require 'fileutils'
require 'net/http'
require 'uri'

BASE = 'https://www.csc.kth.se/cvap/actions/'
SRC = File.join(__dir__, 'data', 'kth')
OUT = File.join(__dir__, 'data', 'curves', 'kth')
WIDTH = 160
HEIGHT = 120
PERSONS = (1..25).freeze
ACTIONS = %w[walking jogging running boxing handwaving handclapping].freeze

# Streams the body of `url` into `target` through a partial file.
def download(url, target)
  partial = "#{target}.part"
  Net::HTTP.get_response(URI(url)) do |response|
    raise "#{url}: HTTP #{response.code}" unless response.is_a?(Net::HTTPSuccess)

    File.open(partial, 'wb') { |io| response.read_body { |chunk| io.write(chunk) } }
  end
  File.rename(partial, target)
end

def fetch(name)
  path = File.join(SRC, 'zip', name)
  return path if File.exist?(path)

  FileUtils.mkdir_p(File.dirname(path))
  puts "downloading #{BASE}#{name}"
  download(BASE + name, path)
  path
end

# First frame range of each sequence: "person01_boxing_d1" => 1..95.
def first_ranges
  File.readlines(fetch('00sequences.txt')).filter_map do |line|
    m = line.match(/^(person\d+_\w+_d\d)\s+frames\s+(\d+)-(\d+)/)
    m && [m[1], m[2].to_i..m[3].to_i]
  end.to_h
end

# Grayscale frames `range` (1-based) of the video, as Float64 bytes in [0, 1].
def frame_bytes(avi, range)
  filter = "select=between(n\\,#{range.first - 1}\\,#{range.last - 1})"
  raw = IO.popen(['ffmpeg', '-v', 'error', '-i', avi, '-vf', filter, '-fps_mode', 'passthrough',
                  '-f', 'rawvideo', '-pix_fmt', 'gray', '-'], 'rb', &:read)
  frames, rest = raw.bytesize.divmod(WIDTH * HEIGHT)
  abort "#{avi}: ffmpeg returned #{raw.bytesize} bytes" unless rest.zero? && frames.positive?
  [raw.unpack('C*').map { |v| v / 255.0 }.pack('E*'), frames]
end

def unpack_videos(action)
  dir = File.join(SRC, 'avi')
  names = PERSONS.map { |p| format('person%02d_%s_d1_uncomp.avi', p, action) }
  missing = names.reject { |n| File.exist?(File.join(dir, n)) }
  return if missing.empty?

  FileUtils.mkdir_p(dir)
  system('unzip', '-o', '-q', fetch("#{action}.zip"), *missing, '-d', dir, exception: true)
end

actions = ARGV.empty? ? %w[walking boxing] : ARGV
bad = actions - ACTIONS
abort "unknown action #{bad.join(', ')}; choose from #{ACTIONS.join(' ')}" unless bad.empty?
system('ffmpeg', '-version', out: File::NULL, exception: true)

ranges = first_ranges
actions.each { |a| unpack_videos(a) }
FileUtils.mkdir_p(OUT)
rows = PERSONS.flat_map do |p|
  actions.filter_map do |action|
    seq = format('person%02d_%s_d1', p, action)
    next warn("#{seq}: no frame range listed; skipped") unless ranges[seq]

    bytes, frames = frame_bytes(File.join(SRC, 'avi', "#{seq}_uncomp.avi"), ranges[seq])
    name = format('p%02d_%s', p, action)
    File.binwrite(File.join(OUT, "#{name}.f64"), bytes)
    [name, action, WIDTH * HEIGHT, frames, "#{name}.f64"]
  end
end
File.open(File.join(OUT, 'index.csv'), 'w') do |io|
  io.puts 'name,label,dims,vertices,file'
  rows.each { |r| io.puts r.join(',') }
end
puts "#{rows.size} curves written to #{OUT}"

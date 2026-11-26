#!/usr/bin/env ruby
# frozen_string_literal: true

# Downloads the UCI "Condition monitoring of hydraulic systems" dataset (Helwig,
# Pignanelli and Schütze, 2015; CC BY 4.0) into example_higher/data/hydraulic/ and
# extracts the six pressure sensors PS1..PS6. Each PSk.txt has 2205 rows (machine
# cycles) of 6000 values (60 s at 100 Hz); hydraulic.jl turns each sensor into a curve
# with 2205 vertices in R^6000, as in Meintrup, Munteanu and Rohde (NeurIPS 2019).
# Files already present are kept.
#
# Usage: ruby example_higher/download_hydraulic.rb

require 'fileutils'
require 'net/http'
require 'uri'

URL = 'https://archive.ics.uci.edu/static/public/447/condition+monitoring+of+hydraulic+systems.zip'
DIR = File.join(__dir__, 'data', 'hydraulic')
ZIP = File.join(DIR, 'hydraulic.zip')
MEMBERS = (1..6).map { |k| "PS#{k}.txt" }.freeze
COLUMNS = 6000

# Streams the body of `url` into `target` through a partial file.
def download(url, target)
  partial = "#{target}.part"
  Net::HTTP.get_response(URI(url)) do |response|
    raise "#{url}: HTTP #{response.code}" unless response.is_a?(Net::HTTPSuccess)

    File.open(partial, 'wb') { |io| response.read_body { |chunk| io.write(chunk) } }
  end
  File.rename(partial, target)
end

# Checks that a sensor file starts with a row of COLUMNS tab-separated numbers.
def check_sensor(path)
  first = File.open(path, &:readline).split("\t")
  return if first.size == COLUMNS && first.all? { |v| Float(v, exception: false) }

  raise "#{path}: expected #{COLUMNS} numbers per row, got #{first.size} fields"
end

FileUtils.mkdir_p(DIR)
missing = MEMBERS.reject { |m| File.exist?(File.join(DIR, m)) }
if missing.any?
  unless File.exist?(ZIP)
    puts "downloading #{URL}"
    download(URL, ZIP)
  end
  system('unzip', '-o', '-q', ZIP, *missing, '-d', DIR, exception: true)
end
MEMBERS.each do |m|
  check_sensor(File.join(DIR, m))
  puts format('%<m>-8s %<mb>6.1f MB', m: m, mb: File.size(File.join(DIR, m)) / 1e6)
end
puts "pressure sensors ready in #{DIR}"

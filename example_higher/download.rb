#!/usr/bin/env ruby
# frozen_string_literal: true

# Downloads the motion-capture trials listed in trials.txt from the CMU Graphics Lab
# Motion Capture Database into example_higher/data/ (ignored by git). Each subject
# contributes its skeleton (.asf) and each trial its motion (.amc). Files already
# present are kept, so an interrupted run can simply be repeated.
#
# Usage: ruby example_higher/download.rb

require 'fileutils'
require 'net/http'
require 'uri'

BASE_URL = 'http://mocap.cs.cmu.edu/subjects'
HERE = __dir__
DATA_DIR = File.join(HERE, 'data')
MANIFEST = File.join(HERE, 'trials.txt')
TRIAL_ID = /\A(\d\d)_(\d\d)\z/
# Text every genuine file of the given kind contains; an HTML error page does not.
SIGNATURE = { '.asf' => ':bonedata', '.amc' => ':DEGREES' }.freeze

def trial_ids
  File.readlines(MANIFEST, chomp: true).filter_map do |line|
    id = line.sub(/#.*/, '').split.first
    next if id.nil?
    raise "#{MANIFEST}: bad trial id #{id.inspect}" unless TRIAL_ID.match?(id)

    id
  end
end

# The files of the given trials, as [subject, file name] pairs, skeletons first.
def wanted_files(ids)
  subjects = ids.map { |id| id[0, 2] }.uniq
  subjects.map { |s| [s, "#{s}.asf"] } + ids.map { |id| [id[0, 2], "#{id}.amc"] }
end

def fetch(url)
  response = Net::HTTP.get_response(URI(url))
  raise "#{url}: HTTP #{response.code}" unless response.is_a?(Net::HTTPSuccess)

  response.body
end

# Downloads one file through a partial file in the data directory, renamed into place
# only after its content is checked.
def download(subject, name)
  target = File.join(DATA_DIR, name)
  return :present if File.exist?(target)

  body = fetch("#{BASE_URL}/#{subject}/#{name}")
  signature = SIGNATURE.fetch(File.extname(name))
  raise "#{name}: not a motion-capture file (no #{signature})" unless body.include?(signature)

  partial = "#{target}.part"
  File.binwrite(partial, body)
  File.rename(partial, target)
  :downloaded
end

FileUtils.mkdir_p(DATA_DIR)
files = wanted_files(trial_ids)
counts = Hash.new(0)
files.each_with_index do |(subject, name), i|
  status = download(subject, name)
  counts[status] += 1
  puts format('[%<i>2d/%<n>d] %<name>-10s %<status>s', i: i + 1, n: files.size, name: name, status: status)
end
puts "#{counts[:downloaded]} downloaded, #{counts[:present]} already present, in #{DATA_DIR}"

#!/usr/bin/env ruby
# frozen_string_literal: true

# Bit-for-bit regression check between a git revision and the working tree.
#
#   ruby test/snapshots/compare.rb BASE_REV
#
# Builds both snapshots (golden.jl, golden_ext.jl) at BASE_REV, in a temporary git
# worktree, and at the current working tree, then reports every entry that differs.
# Exit status is 0 only when both snapshots are identical. The scripts always come
# from the working tree, so BASE_REV must expose the same API they call.
require 'tmpdir'
require 'open3'

ROOT = File.expand_path('../..', __dir__)
SNAP = File.join(ROOT, 'test', 'snapshots')
SCRIPTS = %w[golden golden_ext].freeze

def run(*cmd, chdir: ROOT)
  out, status = Open3.capture2e(*cmd, chdir: chdir)
  abort "#{cmd.join(' ')} failed:\n#{out}" unless status.success?
  out
end

def snapshot(project, script, out)
  run('julia', "--project=#{project}", File.join(SNAP, "#{script}.jl"), out)
end

abort 'usage: compare.rb BASE_REV' unless ARGV.size == 1

ok = true
Dir.mktmpdir('frechet-snap') do |tmp|
  base = File.join(tmp, 'base')
  run('git', 'worktree', 'add', '--detach', base, ARGV[0])
  begin
    run('julia', "--project=#{base}", '-e', 'using Pkg; Pkg.instantiate()')
    SCRIPTS.each do |script|
      snapshot(base, script, File.join(tmp, "#{script}_base.jls"))
      snapshot(ROOT, script, File.join(tmp, "#{script}_new.jls"))
      cmp = ['julia', File.join(SNAP, 'compare.jl'), File.join(tmp, "#{script}_base.jls"), File.join(tmp, "#{script}_new.jls")]
      out, status = Open3.capture2e(*cmp)
      puts "== #{script}", out
      ok &&= status.success?
    end
  ensure
    run('git', 'worktree', 'remove', '--force', base)
  end
end
exit(ok ? 0 : 1)

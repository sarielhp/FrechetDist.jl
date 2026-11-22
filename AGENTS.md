# Agent instructions

## Quality gate

- `tools/gate` (or `tools/check`) is the fast gate: `julia-audit --strict`,
  `julia-static-analysis --syntax --format`, the temporary-file policy, and the test
  suite (`Pkg.test()`, which includes Aqua). `tools/gate --full` adds JET and
  ExplicitImports; run it (or `tools/commit --full`) for changes to `src/`.
- `tools/commit [--full] <message>` runs the gate unless the tree already passed it,
  stages everything not ignored, and commits.
- Numerical changes must be checked bit for bit: `ruby test/snapshots/compare.rb REV`
  builds both snapshots at `REV` and in the working tree and lists differing entries.
  A change that alters results belongs in its own commit and must be flagged.
- Example and benchmark dependencies live in `examples/Project.toml`, never in the
  package's `Project.toml`.

## Releasing

- **"release"**: releases are frequent and do not register the package. Add a changelog
  entry in `docs/src/changelog.md` first, then run `tools/release [--no-install] <message>`
  from `main`: it commits pending work behind the gate, bumps the patch version in
  `Project.toml` behind the full gate, and pushes, verified. For a minor or major
  version, commit with `tools/commit` and run `tools/bump minor` (or `major`).
- **"bump"**: `tools/bump [major|minor|patch]` on a clean tree; `--resume` finishes an
  interrupted bump, `--abandon` drops one that published nothing.
- **"publish"**: registration is rare. Run `tools/publish` only after major changes or
  important breaking changes, and only when explicitly asked; ordinary releases are not
  registered. When asked, run `tools/publish [--dry-run] [release notes]` after a
  release: it waits for the CI workflow on the pushed commit and requests
  JuliaRegistrator registration. The registry PR auto-merges (about 15 minutes for a new
  version) and TagBot creates the tag.
- The `julia` compat bound must stay at a version where the package resolves and its tests
  pass (currently 1.10); AutoMerge installs on the lowest compatible Julia.

## Julia Tooling and Authorization

Follow the Julia [GUIDELINES](~/prog/standards/julia/GUIDELINES.md) and
[RATIONALE](~/prog/standards/julia/RATIONALE.md). Code is formatted with
JuliaFormatter's default style (`.JuliaFormatter.toml`); formatting changes go in their
own commit.

Use the exact command forms when applicable:

| Command | Use |
|---|---|
| `rtk proxy julia-audit --strict` | Depth, complexity, size, idiom and trap audit |
| `rtk proxy julia-static-analysis [paths]` | Syntax, formatting; JET and imports for a project |
| `rtk proxy julia-static-analysis --aqua` | Aqua package checks (also run by the test suite) |
| `rtk proxy julia --project=. -e 'using Pkg; Pkg.test()'` | Run the test suite |

Projects that copied the lifecycle templates also have these (`tools/<name> --help`
says what each one does and runs nothing):

| Command | Use |
|---|---|
| `rtk proxy tools/build` | Check the declared scripts and precompile |
| `rtk proxy tools/gate` or `rtk proxy tools/check` | Fast gate; add `--full` for release verification |
| `rtk proxy tools/lint` / `tools/test` | Static feedback / the test suite |
| `rtk proxy tools/snapshot [message]` | Local checkpoint without gate/push |
| `rtk proxy tools/commit [--full] <message>` | Requested gated local commit |
| `rtk proxy tools/bump [major\|minor\|patch]` | Requested version bump and verified push |
| `rtk proxy tools/bump --resume` | Resume an inspected failed push without bumping again |
| `rtk proxy tools/release [--no-install] <message>` | Requested release: commit, bump, push, install to declared targets unless `--no-install`; never registers |
| `rtk proxy tools/publish [release notes]` | Registration in General; see below |

**Publishing is rare and needs an explicit request.** A release pushes a new version
and is routine once asked for; it does not register the package. `tools/publish`
registers the released version in Julia's General registry, which is public and
cannot be withdrawn. Run it only when the user explicitly asks to publish or
register, typically after major or breaking changes, and never as part of a
release, a bump, or "finishing up". `tools/publish --dry-run` shows what would be
registered and changes nothing.

Run routine audits and tests for authorized work. Existing command approval allows
execution; it does not authorize commits, releases, registration, or changes outside
the task. Reuse authorization already given in the session. Keep sandboxing and manual
approval review enabled. Bootstrap the analysis packages separately, when authorized, with
`julia-install-tools`; the gates never install anything. Never enlarge an allow marker
(`# julia-audit: allow-<rule>`) or add entries to `.julia-audit-baseline` to hide a new
finding; remove baseline entries for code you fix.

Before any refactor of numerical code, record the outputs on fixed inputs and require
bit-identical results afterwards; a result change is its own deliberate commit.

Create every temporary file under `<tmpdir>/<user>/<project>/` (standards README,
"Temporary and cached files"): through `ProjectTmp` (`tools/tmp_space.rb`) in Ruby, or the
language's call with an explicit base under that tree. Never use the system
temporary directory directly or name a literal `/tmp` path; the gate rejects it
unless the line says `tmp-exempt: <reason>`. Files reused across runs are caches,
not temporary files: put them under `$XDG_CACHE_HOME/<project>/`.

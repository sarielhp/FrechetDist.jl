# Agent instructions

## Quality gate

- `tools/check` runs the static analysis (syntax, JuliaFormatter default style, JET,
  ExplicitImports), `julia-audit` (code metrics) and the tests. All must pass.
- `tools/commit <message>` runs the gate, stages everything not ignored, and commits.
- Numerical changes must be checked bit for bit: `ruby test/snapshots/compare.rb REV`
  builds both snapshots at `REV` and in the working tree and lists differing entries.
  A change that alters results belongs in its own commit and must be flagged.

## Releasing

- **"release"** / **"bump"**: releases are frequent and do not register the package. Run
  `tools/release [--minor|--major] [--dry-run] <message>` from `main`: it commits pending
  work behind the gate, bumps the version in `Project.toml` (patch by default, via
  `tools/bump`), pushes, and verifies the remote. Add a changelog entry in
  `docs/src/changelog.md` first.
- **"publish"**: registration is rare. Run `tools/publish` only after major changes or
  important breaking changes, and only when explicitly asked; ordinary releases are not
  registered. When asked, run `tools/publish [--dry-run] <release notes>` after a release
  from a clean, pushed `main`: it requests JuliaRegistrator registration. The registry PR
  auto-merges (about 15 minutes for a new version) and TagBot creates the tag.
- The `julia` compat bound must stay at a version where the package resolves and its tests
  pass (currently 1.10); AutoMerge installs on the lowest compatible Julia.

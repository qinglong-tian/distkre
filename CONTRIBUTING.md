# Contributing

Changes to the maintained library belong in `R/`, with matching help updates
in `man/` and mathematical checks in `tests/`. Keep the method's normalization,
default unclipped predictions, fixed-update semantics, and client/server
information boundaries explicit.

Run package installation, package checks, examples, and the reference workflow
before submitting a change. New estimator changes should be checked against
independent calculations, not only against the implementation itself.

The numerical engine and reference metrics under `reproduction/` are frozen
paper artifacts. Do not silently update them when changing the public API.
Use a new protocol/configuration and retain the old references when numerical
settings change. Document any portability patch in `source-provenance.json`.

Never commit downloaded data, local paths, caches, fitted models, temporary
outputs, or credentials. `reproduction/release-files.sha256` identifies the
publication snapshot; update it intentionally when preparing a new release.

# distkre

Distributed density-ratio estimation with shared random Fourier features.
The package estimates `q(x) / p(x)`, where `q` is the target density and
`p = sum(pi_c * p_c)` is the pooled source density. Source clients may have
different distributions and sample sizes.

There are two entry points:

- **Use the method:** install the R package and follow the examples below.
- **Reproduce the paper:** use the independent workflows in
  [reproduction/](reproduction/README.md). Saved per-replication metrics and
  frozen numerical code are included; raw datasets are downloaded separately.

## Installation

Download or clone this repository, open a terminal in its root, and run:

```sh
R CMD INSTALL .
```

The package requires R >= 4.1 and uses only R's standard libraries. No Python
installation is needed to use the estimator. Reproducing the application
experiments additionally requires Python >= 3.9 for downloading and checking
public data; it uses only Python's standard library.

## Fit a density ratio

```r
library(distkre)

set.seed(42)
sources <- list(matrix(rnorm(200, -0.6), ncol = 2),
                matrix(rnorm(300, 0.6), ncol = 2))
target <- matrix(rnorm(240, 0.2), ncol = 2)

fit <- dkre_rff(sources, target,
                lambda = 0.1, n_features = 64,
                method = "gd", iterations = 1000, seed = 42)
ratio <- predict(fit, target)
```

Use numeric matrices with the same columns, ordering, and preprocessing at
every site. Supply one matrix per source client. The target sample supplies
the numerator distribution. See [the user guide](docs/user-guide.md) for
parameter selection, preprocessing, weights, and saved models.
The [implementation notes](docs/algorithm.md) give the mathematical conventions.

Always specify `method = "gd"` for Algorithm 1. The API also provides
`method = "closed"`, a numerical reference that aggregates second moments.
The default remains `"closed"` for compatibility with earlier versions.

## Separate clients and server

`dkre_client()`, `dkre_client_update()`, and `dkre_server_fit()` expose the
client/server interface. Each client retains its source observations and
second-moment matrix; it returns a vector for each server update.
[examples/client-server.R](examples/client-server.R) shows a complete local
simulation of this interface. A deployment supplies its own transport
callback; the package does not start remote processes.

```sh
Rscript examples/quickstart.R
Rscript examples/client-server.R
Rscript examples/tuning.R
```

## Reproduce the paper

From the repository root:

```sh
# Rebuild and verify the manuscript tables, supporting table, and figures.
# Offline: this uses archived metrics and performs no model fitting.
Rscript reproduction/reproduce.R --mode=reference

# Download public data, refit a small check collection, and render its outputs.
# Application settings are unchanged; simulation sizes are reduced.
Rscript reproduction/reproduce.R --mode=smoke

# Refit the complete production study; allow substantial computation time.
Rscript reproduction/reproduce.R --mode=full --workers=4
```

Outputs go to `reproduction/output/<mode>/`. They are not committed.
[The reproduction guide](reproduction/README.md) explains every paper table,
seed, data source, frozen setting, expected count, and computational limit.
The paper applications are Heart Disease, Beijing Air Quality, and
positive-unlabeled Spambase; exploratory applications are not part of this
public release.

## Repository layout

```text
R/                  maintained estimation library
man/                R help and API examples
examples/           complete usage and tuning examples
tests/              mathematical and client/server checks
docs/               practitioner guide and contribution instructions
reproduction/       reproducible paper workflows and frozen numerical engine
  config/           frozen paper simulation configuration and counts
  reference/        archived per-replication metrics and reference summaries
  engine/           historical numerical computations, isolated from the API
  tests/            replay checks against the archived paper metrics
paper/              numerical-section LaTeX and dataset references
.github/workflows/  automated package and reproduction checks
```

## Checks, citation, and license

```sh
R CMD build .
R CMD check --no-manual distkre_0.1.0.tar.gz
python3 reproduction/verify-release.py
```

Use `citation("distkre")` for the software citation. Metadata is also supplied
in [CITATION.cff](CITATION.cff). Code is licensed under GPL version 3 or later;
datasets retain their original terms and are not bundled. The maintainer is
Qinglong Tian (<qltian2021@gmail.com>).

# Reproducing the paper

Run commands from the repository root. R >= 4.1 and Python >= 3.9 are needed
for the complete workflow; no additional R or Python packages are required.
Original computations used R 4.6.0 on macOS with an M1 Pro. Reproduction scripts
explicitly use Mersenne-Twister / Inversion / Rejection RNG settings, and launch
R subprocesses with `--vanilla` to avoid user startup files. See
[the original session record](reference/original-session.txt). Repeated runs
use explicit seeds. Floating-point results can vary across R versions and
math libraries; the replay checks compare numerical values with tolerances.

## Three modes

```sh
Rscript reproduction/reproduce.R --mode=reference
Rscript reproduction/reproduce.R --mode=smoke
Rscript reproduction/reproduce.R --mode=full --workers=4
```

| Mode | Computation | Purpose |
|---|---|---|
| `reference` | Reaggregate archived metrics, verify numerical summaries, regenerate tables and simulation figures | Offline reconstruction of the reported outputs; no fits |
| `smoke` | Five reduced simulation fixtures plus one original-setting fixture for each application | Test the complete pipeline; simulation outputs are not paper results |
| `full` | Refit all 3,306 production simulation datasets and 100 application fixtures | Recreate the paper computations from generated data and downloaded observations |

Every mode creates `PAPER_TABLES.md` and six table fragments in `tables/`
(five manuscript tables and one supporting unequal-size comparison).
Reference and full modes also create simulation figures under
`simulation/summary/production/`; smoke and subset runs omit publication plots.
Tables contain
tabular environments, ready to insert into the manuscript. Complete helper
summaries and job files are retained alongside them. Reference reconstruction
sets unavailable per-job elapsed times to zero: use the original metrics for
historical timings, and new runs for timings on your machine. Timings are not
part of the numerical replay equality checks.

Full runs are substantial: exact kernel systems at source sizes 4,000 and
8,000 are cubic-cost references, and each candidate has an optimization
budget. Start with smoke mode, choose a worker count for available memory,
and expect a full run to take many hours or longer. Four workers is an example,
not a memory requirement or runtime guarantee. Windows executes sequentially;
macOS and Linux can use multiple workers. No Apple-specific software, BLAS
switching, shell supervisor, or parallel R extension is required.

## Paper-output map

| Paper output | Reconstructed file | Metric input / computation |
|---|---|---|
| Simulation accuracy by source/target size | `tables/simulation_accuracy.tex` | Seven accuracy settings, independent-test RMSE and MCSE |
| Accuracy--communication comparison | `tables/simulation_communication.tex` | First recorded checkpoint within 2% of the RF mean at lambda 0.01 |
| Client heterogeneity | `tables/simulation_heterogeneity.tex` | Pooled truth fixed; three heterogeneity levels |
| Supporting unequal-size comparison (discussed in the text) | `tables/simulation_imbalance.tex` | 70/10/10/10 source proportions; uniform-weight ablation retained |
| Hospital and station prediction | `tables/applications_natural.tex` | Four held-out hospitals and twelve held-out stations |
| Positive-unlabeled spam ranking | `tables/applications_spam.tex` | Twenty splits, random and heterogeneous client ownership |
| Three simulation figures and KDE diagnostic | `simulation/summary/production/study*.png` and `.svg` | Same per-replication inputs and frozen plotting code |

The readable LaTeX sections are in `paper/`. The application tables exclude
Local-J, as in the current manuscript. Other diagnostic rows in the archived
metrics are preserved but not presented as additional positive applications.

## Frozen simulation design

`config/simulation.R` is a readable R object specifying all seeds, grids,
budgets, and final case counts. Its production counts are:

| Case | Replications |
|---|---:|
| N=1000, M=1000 | 500 |
| N=2000, M=1000 | 306 |
| N=4000, M=1000 | 200 |
| N=8000, M=1000 | 200 |
| N=2000, M=250 | 500 |
| N=2000, M=500 | 500 |
| N=2000, M=2000 | 200 |
| Unequal source sizes | 500 |
| Each of three heterogeneity levels | 100 |
| Communication trajectories | 100 |

Final counts were determined from an independent 240-dataset pilot, targeting
paired MCSE 0.0025 with a cap of 500. Production reruns use these already frozen
counts. The full-kernel extension is included at N=4000 and N=8000 on the same
400 production datasets; it does not add independent datasets. The original
extension used Accelerate only for speed. The portable rerun uses the same
deterministic data and validation grids and includes the reference in each job.

To reproduce the pilot separately:

```sh
Rscript reproduction/run-simulations.R --mode=pilot --output=reproduction/output/pilot/simulation
```

To replay one production fixture at its original settings:

```sh
Rscript reproduction/run-simulations.R --mode=full \
  --case=accuracy_N1000_M1000 --replications=1 \
  --output=reproduction/output/single/simulation
```

Subset runs do not recreate the full paper table. Their summaries report the
number completed relative to the frozen count. Use a separate output directory
for each configuration.

## Application datasets and settings

The downloader retrieves only these original UCI datasets:

- [Heart Disease](https://archive.ics.uci.edu/dataset/45/heart+disease): four
  processed hospital files, 920 patients, observed disease labels. Leave each
  hospital out; five seeds (1--5).
- [Beijing Multi-Site Air Quality](https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data):
  twelve original station CSVs. The adapter creates 17,164 station-days,
  requiring at least 18 measurements per daily mean and an observed response.
  All twelve targets, five subsampling seeds (1--5), chronological calendar-month
  roles and one-month gaps.
- [Spambase](https://archive.ics.uci.edu/dataset/94/spambase): 4,601 emails and
  57 features. Seeds 201--220. Positive-only labeling and provider allocation
  are constructed on observed records; these are not observed operational
  provider assignments. Random and heterogeneous ownership use the same records.

`datasets.json` pins SHA-256 hashes of every required raw file. The downloader
validates bytes before replacing any dataset files and records retrieval
provenance. An upstream change is an error, not a silent new benchmark.
Raw data, derived caches, fitted models, and local results are ignored by Git.
Dataset citations are in `paper/application-references.bib`; dataset terms
remain separate from the software license.

```sh
python3 reproduction/fetch-data.py
python3 reproduction/fetch-data.py --check
Rscript reproduction/run-applications.R --mode=full --workers=4
```

For a single family, pass `--datasets=heart`, `--datasets=beijing`, or
`--datasets=spambase` and a new `--output` directory. The default is all three.
Hospital/station ratios use R=128, 3,000 updates, lambda {0.01,0.1}, and three
bandwidth multipliers. Source roles are capped at 80 fitting/80 validation/160
model-training observations; target roles at 240 fitting/80 validation/480
evaluation observations. Spam uses R=512, 5,000 updates, lambda
{0.001,0.01,0.1}, and the same three bandwidth multipliers. Local fitting has
common target access and independent source-only tuning. The local prediction
comparison includes fitting separate predictors; its loss is a weighted average
of their target losses. The spam comparison uses each client's own evaluation
records and reports global and within-provider ranking separately.

Hospital and Beijing shifts retain original domains and outcomes. Conditional
invariance is not established. Weighting does not improve Beijing's aggregate
prediction; the hospital average gain is concentrated in Switzerland. These
findings and the sample caps are preserved in the manuscript and references.

## Checks and provenance

```sh
python3 reproduction/verify-release.py
Rscript reproduction/tests/check-replay.R
```

The replay check refits one full-setting simulation dataset and the three
smoke application fixtures and compares their outputs with the archived metrics.
Reference mode independently checks all summary RMSEs/MCSEs and all application
summary values. Full mode additionally compares every saved production metric
with the archived collection (tolerance 1e-8; timings excluded), and requires
complete job coverage. Smoke mode checks its original-setting application
fixtures against their archived results. To verify a full application-only run:

```sh
Rscript reproduction/verify-results.R --scope=applications --input=reproduction/output/full-app-check
```

The package tests verify the mathematical formulas, finite GD
trajectories, vector-only transport, clipping, stopping, and saved-model behavior.
The frozen engine checks additionally verify data generation and parameter
selection. A changed configuration or source checksum prevents resuming into
an existing output directory. Job files are saved atomically; saved failures
are reported and stop the run rather than being omitted from summaries.

`source-provenance.json` maps original numerical helpers to the release copies.
The core estimator is unchanged. Portable orchestration replaces machine-specific
supervisors. The data adapter changes only the Spambase path, the analysis helper
accepts a portable output directory and can omit plots for incomplete subsets,
and natural jobs use atomic saves. `collect-applications.R` separately audits
and aggregates completed application jobs without refitting them.
`release-files.sha256` covers the files shipped in this release.

## Exploratory cases

The original research workbench is preserved outside this clean repository.
Turbine CO/NOx, policy evaluation, Adult risk auditing, ecological examples,
and other screens remain available there for future work. They are deliberately
outside the paper reproduction commands and are not required to recreate any
reported table in the current draft.

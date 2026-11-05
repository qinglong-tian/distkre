# Publication validation

Validated locally on November 5, 2026, using R 4.6.0 on macOS.

| Check | Result |
|---|---|
| R package build and installation | Passed |
| `R CMD check --no-manual` | Status OK: no errors, warnings, or notes |
| Practitioner quickstart, client/server, and tuning examples | All passed |
| Frozen simulation engine mathematical/RNG/partition checks | Passed |
| Positive-unlabeled evaluation-label isolation and allocation checks | Passed |
| Offline reconstruction | All 40,942 simulation metric rows reaggregated; all archived summary RMSEs and MCSEs verified |
| Application reference reconstruction | All archived summary values verified |
| Fresh production simulation fixture | N=1000, M=1000, replication 1; original settings matched archived numerical metrics within 1e-8 |
| Fresh complete application collection | All 80 hospital/station fixtures and 20 spam fixtures refitted from downloaded original data |
| Application numerical replay | All 1,200 hospital/station metric rows and 960 spam metric rows matched exactly (maximum absolute difference 0) |
| Application domain, row, and temporal-block isolation | Passed for all saved fixtures |
| Complete smoke entry point | Passed, including reduced simulations, three original-setting applications, rendering, and numerical verification |

The complete 3,306-dataset production simulation study was not refitted during
this organization work. Its full rerun is scripted; local validation used the
entire archived collection, one fresh production dataset, and reduced simulation
fixtures covering all study types. Reference mode reconstructs results from
archived metrics and is explicitly separate from a fresh fit.

GitHub Actions is configured for package checks on Linux, macOS, and Windows,
and for reference reconstruction and fixture replay on Linux. Those hosted
jobs will run when the repository is uploaded; local testing does not assert
that the hosted jobs have already executed.

The existing GPL version 3-or-later licensing is retained. Maintainer metadata
uses Qinglong Tian and qltian2021@gmail.com, as supplied by the author. Raw
datasets and generated local files are excluded from the publication snapshot.
The original workbench and all unused exploratory cases remain preserved
outside this repository.

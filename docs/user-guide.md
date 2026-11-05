# Using distkre

## Inputs and estimand

Each source matrix has observations in rows and covariates in columns.
`target_data` has the same columns. The estimator targets

```text
r(x) = q(x) / sum_c pi_c p_c(x),   pi_c = n_c / sum_d n_d.
```

One source is supported as a special case. For heterogeneous clients, fitting
separately at each client estimates different ratios `q/p_c`; averaging those
functions generally does not recover the pooled ratio.

Prepare numeric features before calling the package. It rejects missing or
non-finite values. Encode categories consistently and apply the same centering
and scaling at all sites. Choose preprocessing using fitting data only. In the
paper applications, target fitting means and standard deviations are shared;
held-out outcomes are never used. Other deployments may use agreed public
scales or a suitable distributed preprocessing procedure.

Target coverage by the pooled source distribution is needed for importance
weighting. A fitted finite value cannot establish coverage where source
observations are absent. Natural shifts can also change outcomes conditional
on the measured covariates, which covariate reweighting alone does not correct.

## Parameters

| Argument | Meaning | How to choose |
|---|---|---|
| `lambda` | Positive regularization on the paper's scale | Compare candidates on held-out ratio loss |
| `bandwidth` | Gaussian kernel length scale | Start with a target-distance scale; validate alternatives |
| `n_features` | Shared random-feature dimension | Balance approximation accuracy against memory and vector size |
| `iterations` | GD update budget | Compare against the closed-form reference during development or inspect diagnostics |
| `seed` | Shared feature randomness | Fix for reproducibility; clients must use the same specification |
| `tolerance` | Optional relative-gradient stopping threshold | Zero reproduces the fixed-update algorithm |

`examples/tuning.R` demonstrates the held-out criterion

```text
0.5 * sum_c pi_c * mean(r(source_validation_c)^2)
    - mean(r(target_validation)).
```

This criterion needs covariates, not target outcome labels. Splits should respect
subjects, repeated measurements, and time dependence where applicable. The
example uses independent simulated data; the paper scripts specify separate
clinical and temporal splits.

## Predictions and downstream weights

`predict(fit, new_matrix)` returns density-ratio estimates, not class
probabilities. Unconstrained KuLSIF can return negative values. Default fits
and predictions remain untruncated to match the estimation method. For
applications requiring nonnegative importance weights, use
`pmax(0, predict(fit, new_matrix))` and document any normalization or clipping.
The paper normalizes downstream weighted losses by total weight. Do not divide
by zero if a fitted sample has no positive weight mass.

The optional `clip = TRUE` truncates negative predictions. Supplying
`upper_bound = B` additionally clips to `[0, B]`. These choices alter the
estimator and are not the package default. The paper's spam ranking uses raw
scores; hospital and station prediction truncate weights below zero.

## Client/server interface

Create one `rff_features()` specification and distribute it to all clients.
Each client creates a `dkre_client()` state. The server's callback broadcasts
the current coefficient vector and returns one `dkre_client_update()` vector
per client, in the same order as `source_sizes`. The server never accepts
client matrices through this callback. See `?dkre_server_fit` and
`examples/client-server.R`.

The client state contains a feature second moment, so memory at a client is
quadratic in the number of features. Messages in each update are linear in
that dimension. Total communication also depends on the number of updates;
vector messages alone do not guarantee a smaller cumulative payload than
one-time matrix aggregation. The local callback example does not measure
network latency or implement an encryption or privacy protocol.

## Diagnostics and saving

`fit$metadata` includes the feature specification, coefficients, iteration
count, stopping reason, and analytical communication count. With fixed updates,
diagnostics refer to the last evaluated iterate, before the final update;
`diagnostics_iteration` identifies it. With early stopping, a final gradient
check may add an exchange.

Save a model with `saveRDS(fit, "ratio-model.rds")` and reload it with
`readRDS()`. Models do not retain raw observations by default; `retain_data = TRUE`
is an explicit opt-in for the convenience wrapper. RDS files should only be
loaded from trusted sources.

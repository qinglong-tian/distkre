# Methodology implementation notes

The maintained estimator is in `R/`. The frozen paper computations are
isolated in `reproduction/engine/`; see `reproduction/README.md` for the
precise settings and result provenance.

## Estimator and conventions

The target is `q(x) / sum_c pi_c p_c(x)`, where `pi_c = n_c / n`.
For a Gaussian kernel `exp(-||x-y||^2 / (2 sigma^2))`, the shared map is

    phi_R(x) = sqrt(2/R) (cos(W_1' x + phase_1), ..., cos(W_R' x + phase_R))

with independent `W_r ~ N(0, sigma^-2 I)` and `phase_r ~ Uniform(0, 2 pi)`.
All sites use the same realizations. Bandwidth, feature count, input dimension,
feature law and random-generation conventions must agree, not merely a seed.

Each client keeps `A_c = Phi_c' Phi_c / n_c` locally. The target server computes
`a = colMeans(Phi_target)`. Equation (13) is

    (sum_c pi_c A_c + lambda I) b = a.

`method = "closed"` solves that system, intentionally aggregating matrices.
The parameter `lambda` is exactly the paper-scale regularization parameter.

## Algorithm 1 and site boundaries

- `dkre_client()` creates the local count and second moment and discards raw
  observations from that state. The state stays at its source site.
- `dkre_client_update()` accepts the broadcast coefficient vector and returns
  only `A_c b`.
- `dkre_server_fit()` accepts target observations, public counts, a shared feature
  specification and a transport callback. The callback returns one R-vector per
  client. The server neither accepts nor aggregates full client matrices.
- `dkre_rff(method = "gd")` simulates these interfaces locally. A real deployment
  supplies the transport callback; network execution is not bundled.

Initialize `b = 0`. At every update the server uses

    g = sum_c pi_c (A_c b) - a + lambda b
    b_next = b - g / (2 + lambda).

The bound `||phi_R(x)||^2 <= 2` justifies the conservative step size.
The default `tolerance = 0` executes exactly T updates and T exchanges, including
when a gradient is zero. Each exchange accounts for `2 C R` scalar coordinates.
With positive tolerance, the relative gradient criterion is

    ||g|| / max(1, ||b|| + ||a||) <= tolerance.

It is checked before updating the current iterate. An early stop therefore
uses one exchange to verify the stopping condition after the preceding update.
The returned metadata reports actual updates and exchanges separately.

Objective and gradient diagnostics use the received vector products. Their
iterate index is explicit; at a fixed-T exit they describe b_(T-1), while the
returned coefficients are b_T. Computing diagnostics this way avoids any extra
matrix transmission or final diagnostic exchange. `trace = TRUE` retains the
history; otherwise only scalar diagnostics are retained.

## Outputs and optional clipping

The default output is `b' phi_R(x)` with no clipping. The optional `clip = TRUE`
applies `max(0,r(x))`; adding `upper_bound = B` applies `min(B,max(0,r(x)))`.
The same rule applies to stored fitted values and `predict()` results.
Kernel KuLSIF also defaults to no clipping and supports the same options.

Random-feature fits keep input dimension and sample counts for prediction and
printing, and omit raw observations by default. The local wrapper has an
explicit `retain_data = TRUE` opt-in. The server-only interface never returns
raw observations or client matrices. These interfaces preserve site boundaries;
they do not assert differential privacy or prevent all inference from updates.

## Bandwidth and parameter selection

For the Gaussian kernel K(x,y) = exp(-||x-y||^2/(2 sigma^2)), our implementation
uses the median pairwise Euclidean distance among target observations as a
bandwidth heuristic when sigma is not supplied. The server shares this bandwidth
with the feature-generation information. A supplied shared feature specification
instead determines the bandwidth. A positive bandwidth must be supplied when
the target median is zero or cannot be computed.

The regularization parameter lambda and number of features R are supplied by
the user. The implementation does not currently select them by cross-validation.
Algorithm 1 uses a fixed number T of updates by default. An optional relative
gradient stopping criterion is available for computational convenience and is
reported separately from the fixed-T procedure. The estimator is unclipped by
default; optional projection onto [0,B] can be requested when B is supplied.

## Kernel reference

`fit_kulsif()` implements equation (6):

    (K_XX/n + lambda I) alpha = -K_XZ 1_m / (n m lambda)
    r(x) = K(x,X) alpha + K(x,Z) 1_m / (m lambda).

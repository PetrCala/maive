# WAIVE: More Aggressive Correction for P-Hacking and Spurious Precision

WAIVE (Weighted Adjusted Instrumental Variable Estimator) provides a
more aggressive correction for p-hacking and spurious precision by
extending MAIVE with exponential-decay weights that downweight both
spuriously precise estimates and extreme outliers.

## Usage

``` r
waive(
  dat,
  method,
  weight,
  instrument,
  studylevel,
  SE,
  AR,
  first_stage = 1L,
  estimate = NULL,
  se = NULL,
  n = NULL,
  study_id = NULL,
  seed = 123
)
```

## Arguments

- dat:

  Data frame with columns bs, sebs, Ns, study_id (optional). Column
  names can be remapped with `estimate`, `se`, `n`, and `study_id`.

- method:

  1 FAT-PET, 2 PEESE, 3 PET-PEESE, 4 EK.

- weight:

  0 no weights, 1 standard weights, 2 MAIVE adjusted weights, 3 study
  weights.

- instrument:

  1 yes, 0 no.

- studylevel:

  Correlation at study level: 0 none, 1 fixed effects, 2 cluster.

- SE:

  SE estimator: 0 CR0 (Huber-White), 1 CR1 (Standard empirical
  correction), 2 CR2 (Bias-reduced estimator), 3 wild bootstrap.

- AR:

  Anderson Rubin corrected CI for weak instruments (available for
  unweighted and MAIVE-adjusted weight versions of PET, PEESE,
  PET-PEESE, not available for fixed effects): 0 no, 1 yes.

- first_stage:

  First-stage specification for the variance model: 1 log (default), 0
  levels. The log stage regresses log(sebs^2) on log(Ns) with a smearing
  retransformation, so the fitted variance is always positive. The
  levels stage regresses sebs^2 on 1/Ns and is the specification in the
  published paper; pass `first_stage = 0` to reproduce it.

- estimate:

  Optional column name to use instead of 'bs'

- se:

  Optional column name to use instead of 'sebs'

- n:

  Optional column name to use instead of 'Ns'

- study_id:

  Optional column name for study identifiers. When not supplied, a
  column named `study_id` is used if present; otherwise, if `dat` has
  four or more columns, the fourth column is used as the study
  identifier and a warning names the column. Columns already mapped
  through `estimate`, `se`, or `n` are never picked: when the fourth
  column is one of them, the single remaining column is used instead,
  and when several remain an error asks for `study_id` (at
  `studylevel = 0` the fallback is skipped silently, since no identifier
  is needed). Any fourth column (a moderator, a year) would otherwise
  drive the study dummies and clustering at every `studylevel` other
  than 0, so name the column explicitly or drop it. When `studylevel` is
  not 0 and no identifier can be found, a warning says the study-level
  options have no effect.

- seed:

  Seed for the wild bootstrap when SE = 3. Use NULL to avoid setting a
  seed (results depend on the current RNG state). Default is 123 for
  historical reproducibility.

## Value

List with the same structure as
[`maive()`](https://petrcala.github.io/MAIVE/reference/maive.md). See
[`?maive`](https://petrcala.github.io/MAIVE/reference/maive.md) for
details.

## Details

Guided, interactive workflow available at https://www.easymeta.org.

For technical details on WAIVE methodology, see:
<https://meta-analysis.cz/waive_ottawa.pdf>

WAIVE combines variance instrumentation (as in MAIVE) with robust
downweighting based on first-stage residuals. Studies with negative
residuals (spurious precision) or extreme residuals (outliers) receive
exponentially reduced influence in the meta-analytic estimate. This
makes WAIVE more aggressive than standard MAIVE at correcting for
p-hacking and handling outliers.

## Examples

``` r
dat <- data.frame(
  bs = c(0.5, 0.45, 0.55, 0.6),
  sebs = c(0.25, 0.2, 0.22, 0.27),
  Ns = c(50, 80, 65, 90)
)

result <- waive(dat,
  method = 3, weight = 0, instrument = 1,
  studylevel = 0, SE = 0, AR = 0
)
#> Warning: Sample size (4) is small for IV estimation. Results may be unreliable. Consider
#> using instrument=0 for small samples.
#> Warning: Very weak instrument detected (F-test = 0.003). Results may be unreliable.
#> Consider using instrument=0 or checking data quality.
```

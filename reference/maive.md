# R code for MAIVE

R package for MAIVE: "Spurious Precision in Meta-Analysis of
Observational Research" by Zuzana Irsova, Pedro Bom, Tomas Havranek,
Petr Cala, and Heiko Rachinger.

## Usage

``` r
maive(
  dat,
  method,
  weight,
  instrument,
  studylevel,
  SE,
  AR,
  first_stage = 0L,
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

  First-stage specification for the variance model: 0 levels, 1 log.

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
  identifier and a warning names the column. Any fourth column (a
  moderator, a year) would otherwise drive the study dummies and
  clustering at every `studylevel` other than 0, so name the column
  explicitly or drop it.

- seed:

  Seed for the wild bootstrap when SE = 3. Use NULL to avoid setting a
  seed (results depend on the current RNG state). Default is 123 for
  historical reproducibility.

## Value

- beta: MAIVE meta-estimate

- SE: MAIVE standard error

- F-test: heteroskedastic robust F-test of the first step instrumented
  SEs

- beta_standard: point estimate from the conventional (non-instrumented,
  inverse-variance weighted) fit of the method chosen

- SE_standard: standard error from the same conventional fit as
  beta_standard

- Hausman: Hausman type test: comparison between MAIVE and standard
  version

- Chi2: 5

- SE_instrumented: instrumented standard errors, one per input row; NA
  for estimates excluded because their fitted variance was not positive

- AR_CI: Anderson-Rubin confidence interval for weak instruments

- pub bias p-value: p-value of test for publication bias / p-hacking
  based on instrumented FAT

- egger_coef: Egger Coefficient (PET estimate)

- egger_se: Egger Standard Error (PET standard error)

- egger_boot_ci: Confidence interval for the Egger coefficient using the
  selected resampling scheme

- egger_ar_ci: Anderson-Rubin confidence interval for the Egger
  coefficient (when available)

- is_quadratic_fit: Details on quadratic selection and slope behaviour

- ek_structure: Structure of the fitted EK model when method=4: "kink",
  "linear", or "intercept" (intercept-only degenerate fit); NA for other
  methods

- boot_result: Boot result

- slope_coef: Slope coefficient

- petpeese_selected: Which model (PET or PEESE) was selected when
  method=3 (NA otherwise)

- peese_se2_coef: Coefficient on SE^2 when PEESE is the final model (NA
  otherwise)

- peese_se2_se: Standard error of the PEESE SE^2 coefficient (NA
  otherwise)

- weights: second-stage weights, one per input row; NA for excluded
  estimates

- instrument_strength: "strong", "weak", "very_weak", "unknown", or
  "not_applicable", from the first-stage F-test

- n_excluded: number of estimates excluded because the first stage
  fitted a non-positive variance for them (0 with the log first stage)

- excluded_rows: positions of the excluded estimates in the input data
  (after completely empty rows are dropped); integer(0) when none

## Details

Guided, interactive workflow available at https://www.easymeta.org.

Data `dat` can be imported from an Excel file via:
`dat <- read_excel("inputdata.xlsx")` or from a csv file via:
`dat <- read.csv("inputdata.csv")` It should contain:

- Estimates: bs

- Standard errors: sebs

- Number of observations: Ns

- Optional: study_id

Default option for MAIVE: MAIVE-PET-PEESE, unweighted, instrumented,
cluster SE, wild bootstrap, AR.

The levels first stage (`first_stage = 0`) regresses the squared
standard errors on 1/N by ordinary least squares, so an individual
fitted variance can be negative even when the intercept is not. Such an
estimate has no instrumented standard error and is excluded from every
quantity built on it: the MAIVE-adjusted weights, the second-stage
regressions (PET, PEESE, PET-PEESE, EK), the first-stage F-test, the
Hausman comparison, the Anderson-Rubin intervals, and the wild
bootstrap. The first stage itself is fitted on all rows; every reported
statistic then uses the same remaining rows. A warning reports the
number of excluded estimates, which is also returned as `n_excluded`
with their positions in `excluded_rows`. The log first stage
(`first_stage = 1`) cannot fit a negative variance, so it never excludes
an estimate.

## Examples

``` r
dat <- data.frame(
  bs = c(0.5, 0.45, 0.55, 0.6),
  sebs = c(0.25, 0.2, 0.22, 0.27),
  Ns = c(50, 80, 65, 90)
)

result <- maive(dat,
  method = 3, weight = 0, instrument = 1,
  studylevel = 0, SE = 0, AR = 0, first_stage = 0
)
#> Warning: Sample size (4) is small for IV estimation. Results may be unreliable. Consider
#> using instrument=0 for small samples.
#> Registered S3 method overwritten by 'clubSandwich':
#>   method    from    
#>   bread.mlm sandwich
#> Warning: Very weak instrument detected (F-test = 0.002). Results may be unreliable.
#> Consider using instrument=0 or checking data quality.
```

# Pre-registration: Meridian GeoX on real Iowa county sales

Written before any analysis was run. Committed to git so the timestamp is checkable.
Every hypothesis below has a test and a pass/fail rule fixed in advance. Results that
come out differently will be reported as they are.

## Source claims (Ye, Sinn, Li, Jiao, Google, 2026)
- C1. Stratified sampling lowers MDE vs random sampling, for every estimator (Sec 5.1.1, Fig 5).
- C2. Stratified sampling + TBR is the best methodology overall (Fig 6). SC/SDID do better when data is noisy (Fig 8).
- C3. Design-aware placebo keeps the false positive rate at the nominal level per study; t-tests inflate it (over 10x at 5%); random placebo is right on median but has a wide spread across studies (Sec 5.2.1, Fig 11).
- C4. Design-aware placebo cuts mean MDE by 26.8% and median MDE by 31.3% vs random placebo (Sec 5.2.3).
- C5. GeoX MDEs are conservative ("simulation-calibrated"), unlike tools that report the best of many candidates (Sec 5.1.2).

## Data
Iowa Liquor Sales (state open data, mirrored on HuggingFace: lmmcfarland1/iowa_liquor_data_csv).
Aggregated to weekly sales in dollars by county (99 counties). These are retailer purchases from
the state wholesaler, not consumer sales. Real geos, real noise, real holidays, real COVID shock.

A "study" is one rolling window of 50 weeks: 26 train, 8 validation, 8 evaluation, 8 test.
Windows start every 8 weeks. Counties with missing weeks in a window are dropped from that window.

## Methods (my implementation, simplified where noted)
- Treatment cell: about 30% of total sales (share between 25% and 35%), as in the paper's single-cell setup.
- Random sampling: random county sets meeting the share constraint.
- Stratified sampling: k-means (k=4) on standardized mean, CV, slope, lag-1 autocorrelation and DTW distance to the average series (paper Sec 3.1), then the same treated fraction drawn inside each stratum, plus the share constraint. Simplification: no Sobol sequence or minimal-discrepancy labels.
- Design search: 2,000 candidates per rule; TBR fit on train, R2 on validation; keep candidates with R2 >= 0.8 (if fewer than 20 pass, keep the top 20); refit on train+validation; drop designs failing an A/A check on evaluation (|z| > 1.645); choose the lowest MDE.
- Estimators: TBR (beta >= 0), synthetic control (nonnegative weights summing to 1), SDID (unit weights with intercept and ridge, nonnegative time weights, Arkhangelsky et al. 2021).
- Effect: cumulative effect per treated county over the 8 test weeks. Model fit on train+validation.
- Studentized statistic: effect / RMSE on evaluation weeks.
- Design-aware placebo: placebo-treated sets drawn from control counties only, with the same rule (same strata, same share constraint within the control pool), filtered to validation R2 >= 0.6, top 200 by R2 kept.
- Random placebo: uniform random sets of the same size from the control pool, no filtering. 200 draws.
- t-test: residual-based z with SE = sigma_residual * sqrt(8), assuming independent weeks.
- MDE = SE * (z_0.95 + z_0.80) / baseline, with SE from the placebo distribution (paper eq. 13), alpha = 0.10 two-sided.

## Hypotheses and decision rules
H1 (C1). Stratified < random on MDE for TBR.
  Test: paired per study, one-sided sign test. PASS if stratified wins in more than half the studies with p < 0.05.

H2 (C2). TBR has the lowest median MDE of the three estimators under stratified sampling and is the best estimator in more studies than either alternative. PASS if both hold.
H2b. SC's advantage over TBR grows with noise. PASS if Spearman correlation between a study's median CV and (MDE_TBR - MDE_SC)/MDE_TBR is positive with p < 0.05.

H3 (C3). A/A false positive rate at nominal 10%, 100 A/A designs per study.
  PASS if: design-aware placebo pooled FPR is within [7%, 13%]; t-test pooled FPR is above 20%;
  interquartile range of per-study FPR is wider for random placebo than for design-aware placebo.

H4 (C4). Design-aware placebo MDE < random placebo MDE on the same designs.
  PASS if median reduction > 0 with sign-test p < 0.05. The size is compared to the paper's 31.3% median, not used for pass/fail.

H5. With a real 5% lift injected into the treated counties' test weeks, the stratified TBR estimate is unbiased and its 90% interval covers.
  PASS if mean relative error is within +/-1 percentage point of the true lift and coverage is within [85%, 95%].

H6 (C5, the attempt to break it). Inject a lift exactly equal to the MDE reported at design time.
  If the MDE is honest, it is detected about 80% of the time.
  PASS if detection rate >= 75%. I expect this one could fail, because the chosen design is the minimum over many candidates, which tends to flatter it.

## Exploratory (no pass/fail)
- FPR in studies whose test window covers March to May 2020 (COVID stockpiling).
- Where Polk County (Des Moines, the largest county) lands in the chosen designs, and whether SC/SDID do better when it is treated.

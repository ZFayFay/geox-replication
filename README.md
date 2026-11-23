# Google's Meridian GeoX, tested on ten years of real Iowa county sales

An independent rebuild of Ye, Sinn, Li and Jiao (Google, 2026), *Meridian GeoX: An Open-Source Framework for Precision and Efficiency in Geo Experiments* ([paper](https://research.google/pubs/pub1090193/), [code](https://github.com/google/meridian-geox)).

The hypotheses and pass/fail rules were committed before any analysis. See [PREREGISTRATION.md](PREREGISTRATION.md) and the first commit. Steps added after seeing a result are labeled exploratory in the notebook.

**Data.** Iowa Liquor Sales (state open data, retailer purchases from the state wholesaler), summed by county. Real geos and real noise; experiments are A/A tests or known lifts injected into the treated counties. 63 rolling 50-week studies, 6,300 A/A tests.

## Results against the pre-registration

| | Hypothesis (paper's claim) | Result | Verdict |
|---|---|---|---|
| H1 | Stratified sampling beats random on MDE (sign test p < 0.05) | 38 of 63 studies, p = 0.065, median MDE 12% lower | Fail (narrowly) |
| H2 | TBR is the most efficient estimator | Median MDE 1.7% vs SDID 2.5%, SC 2.7%; best in 43 of 63 | Pass |
| H2b | SC gains on TBR as noise grows | Spearman -0.21 (wrong sign, p = 0.11) | Fail |
| H3 | Design-aware placebo keeps A/A false positives near 10% | 22% (95% CI 18% to 25%); random placebo 12%; t-test 26% | **Fail** |
| H4 | Design-aware placebo lowers MDE vs random placebo | 42% lower (paper 31%), in all 63 studies | Pass, but see H3 |
| H5 | 5% lift: unbiased, 90% interval covers 85% to 95% | Bias -0.2 points; coverage 44% | **Fail** |
| H6 | Injecting the promised MDE is detected at least 75% of the time | 33% for the chosen design (promised 1.6%, measured 3.7%) | **Fail** |

**Cross-check with Google's library.** `libcheck/geox_lib_aa.py` runs `meridian-geox` 1.0.1 with its defaults on the same data at daily grain: 61 windows, 60 analyzed. It flagged 42% of A/A tests at the 10% level, and its 90% interval covered a 5% lift 58% of the time, while its average estimate (5.3%) was right.

![False positives](figures/fig1_false_positives.png)
![MDE promise](figures/fig3_mde_promise.png)

## Caveats

One public dataset. County liquor orders are lumpy, especially by day, and counties are not media markets. My design search simplifies GeoX's (no Sobol sequences or minimal-discrepancy strata). Google's own evaluation used 192 advertisers' data. The point is not that GeoX is broken, but that its intervals were too narrow on this data, both in my implementation and in theirs, so a backtest on your own history is worth doing.

## Files

- `geox_replication.ipynb`: the full analysis, executed (about 3 minutes on 2 CPU cores).
- `geox.py`: the same code as a script (percent format).
- `iowa_weekly_county.parquet`, `iowa_daily_county.parquet`: the aggregated data.
- `libcheck/`: the run of Google's library and its saved output.
- `figures/results.json`: every number quoted above.

## Run it

```
pip install numpy pandas scipy scikit-learn matplotlib pyarrow
jupyter nbconvert --to notebook --execute geox_replication.ipynb
# optional, about 25 minutes:
pip install meridian-geox && python libcheck/geox_lib_aa.py libcheck/geox_lib_aa_results.json 60
```

The raw data is the Iowa Liquor Sales dataset (data.iowa.gov), mirrored on Hugging Face as `lmmcfarland1/iowa_liquor_data_csv` (7 GB). `libcheck/agg_iowa.py` streams it into the daily county file.

## Credit

GeoX belongs to Steven Ye, Fabian Sinn, Yunxiao Li, Yang Jiao and the Google Meridian team. This is an independent rebuild for learning. Any mistakes are mine.

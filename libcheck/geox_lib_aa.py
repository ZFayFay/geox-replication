"""A/A check of the real Meridian GeoX library (v1.0.1) on Iowa county daily sales.
For each window: run_design on the pre-period, take the top design, analyze the next 28 days with no lift (A/A)
and with a 5% lift injected into the treated counties. Records p-values and intervals."""
import sys, json, time, datetime, warnings, logging
import numpy as np, pandas as pd
import meridian_geox as geox
warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
D = pd.read_parquet("iowa_daily_county.parquet")
D["location"] = D.county.astype(int).astype(str)
full = D.pivot_table(index="date", columns="location", values="sales", aggfunc="sum")
full = full.reindex(pd.date_range(full.index.min(), full.index.max(), freq="D")).fillna(0.0)
PRE, TEST = 168, 28
GAP = int(sys.argv[3]) if len(sys.argv) > 3 else 0   # days of pre-period withheld from the design step
starts = list(range(0, len(full) - PRE - TEST, int(sys.argv[2]) if len(sys.argv) > 2 else 120))
out_path = sys.argv[1]
res = []
for i, s in enumerate(starts):
    t0 = time.time()
    win = full.iloc[s:s + PRE + TEST]
    keep = win.columns[(win.iloc[:PRE] > 0).mean() > 0.2]          # counties with orders on at least 20% of days
    win = win[keep]
    long = win.stack().rename("conversions").reset_index().rename(columns={"level_0": "date"})
    long["spend"] = long.conversions * 0.01
    pre = long[long.date < win.index[PRE - GAP]]
    try:
        cfg = geox.DesignConfig(experiment_duration=datetime.timedelta(days=TEST), experiment_types=geox.ExperimentType.HOLDBACK,
                                methodology=geox.Methodology.TBR, geo_assignment_rule=geox.GeoAssignmentRule.STRATIFIED_SAMPLING,
                                cell_count=1, design_output_count=1, n_candidates=20_000, seed=i)
        ds = geox.run_design(pre, cfg, geox.Constraints(max_conversions_percent=0.3))
        design = list(ds.designs.values())[0]
        treated = sorted(list(design.designs.values())[0].treatment_geos)
        mde = list(design.designs.values())[0].minimum_detectable_effect
        row = dict(window=i, gap=GAP, start=str(win.index[0].date()), test_start=str(win.index[PRE].date()), n_geos=len(keep), n_treated=len(treated), mde=float(mde))
        for lift in [0.0, 0.05]:
            ana = long.copy()
            m = (ana.date >= win.index[PRE]) & ana.location.isin(treated)
            ana.loc[m, "conversions"] *= (1 + lift)
            acfg = geox.AnalysisConfig(design=design, analysis_start_date=win.index[PRE], analysis_end_date=win.index[-1], n_placebo_candidates=20_000)
            r = geox.analyze(ana, acfg)
            m_ = r.results["cell_1"]
            for k in ["lift", "percent_lift"]:
                e = getattr(m_, k)
                row[f"{k}_{lift}"] = dict(point=float(e.point_estimate), lo=float(e.lower_bound), hi=float(e.upper_bound), sd=float(e.standard_deviation), p=float(e.p_value))
        row["sec"] = round(time.time() - t0, 1)
    except Exception as e:
        row = dict(window=i, error=repr(e)[:500])
    res.append(row); json.dump(res, open(out_path, "w"), default=str)
    print(i, row.get("sec"), row.get("error", ""), flush=True)

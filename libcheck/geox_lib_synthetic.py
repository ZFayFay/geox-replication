"""Negative control for the library check: the same GeoX run on synthetic daily data where every geo is exchangeable
by construction (common shocks + independent AR(1) geo noise, sizes as varied as real counties). If GeoX is calibrated here,
the inflated false-positive rate on Iowa comes from the real data, not from how the library is called.
"""

import sys, json, time, datetime, warnings, logging
import numpy as np, pandas as pd
import meridian_geox as geox

warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)
PRE, TEST, G, N_WIN = 168, 28, 50, int(sys.argv[2]) if len(sys.argv) > 2 else 40
GAP = (
    int(sys.argv[3]) if len(sys.argv) > 3 else 0
)  # days of pre-period withheld from the design step
ABS = (
    len(sys.argv) > 4 and sys.argv[4] == "abs"
)  # GeoX option: compare |t| with absolute placebo t-statistics
out_path = sys.argv[1]
res = []
for w in range(N_WIN):
    rng = np.random.default_rng(500 + w)
    T = PRE + TEST
    size = np.exp(rng.normal(9, 1.0, G))  # heterogeneous geo sizes
    t = np.arange(T)
    common = 0.15 * np.sin(2 * np.pi * t / 7) + np.cumsum(
        rng.normal(0, 0.01, T)
    )  # weekly cycle + random walk, shared
    e = np.zeros((G, T))
    eps = rng.normal(0, 0.10, (G, T))
    for k in range(1, T):
        e[:, k] = 0.5 * e[:, k - 1] + eps[:, k]  # independent AR(1) geo noise
    Y = size[:, None] * np.exp(common[None, :] + e)
    dates = pd.date_range("2020-01-01", periods=T, freq="D")
    long = pd.DataFrame(
        {
            "date": np.tile(dates, G),
            "location": np.repeat([str(g) for g in range(G)], T),
            "conversions": Y.ravel(),
        }
    )
    long["spend"] = long.conversions * 0.01
    t0 = time.time()
    try:
        cfg = geox.DesignConfig(
            experiment_duration=datetime.timedelta(days=TEST),
            experiment_types=geox.ExperimentType.HOLDBACK,
            methodology=geox.Methodology.TBR,
            geo_assignment_rule=geox.GeoAssignmentRule.STRATIFIED_SAMPLING,
            cell_count=1,
            design_output_count=1,
            n_candidates=20_000,
            seed=w,
        )
        ds = geox.run_design(
            long[long.date < dates[PRE - GAP]], cfg, geox.Constraints(max_conversions_percent=0.3)
        )
        design = list(ds.designs.values())[0]
        treated = sorted(list(design.designs.values())[0].treatment_geos)
        row = dict(
            window=w,
            gap=GAP,
            abs_stats=ABS,
            n_treated=len(treated),
            mde=float(list(design.designs.values())[0].minimum_detectable_effect),
        )
        for lift in [0.0, 0.05]:
            ana = long.copy()
            m = (ana.date >= dates[PRE]) & ana.location.isin(treated)
            ana.loc[m, "conversions"] *= 1 + lift
            r = geox.analyze(
                ana,
                geox.AnalysisConfig(
                    design=design,
                    analysis_start_date=dates[PRE],
                    analysis_end_date=dates[-1],
                    n_placebo_candidates=20_000,
                    use_absolute_placebo_statistics=ABS,
                ),
            )
            e_ = r.results["cell_1"].percent_lift
            row[f"percent_lift_{lift}"] = dict(
                point=float(e_.point_estimate),
                lo=float(e_.lower_bound),
                hi=float(e_.upper_bound),
                p=float(e_.p_value),
            )
        row["sec"] = round(time.time() - t0, 1)
    except Exception as ex:
        row = dict(window=w, error=repr(ex)[:300])
    res.append(row)
    json.dump(res, open(out_path, "w"))
    print(w, row.get("sec"), row.get("error", ""), flush=True)

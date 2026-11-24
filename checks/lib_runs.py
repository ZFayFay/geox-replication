"""Reruns of Google's meridian-geox library, with warnings recorded.

Same data preparation, windows, seeds and config as libcheck/geox_lib_aa.py and geox_lib_synthetic.py,
with four additions:
  1. library warnings are NOT silenced: every logging.WARNING message is recorded per window;
  2. the number of placebos actually used, and how many of them are distinct, is recorded;
  3. a third analysis injects a lift equal to the MDE the library promised for its chosen design;
  4. n_candidates / n_placebo_candidates / library version (PyPI 1.0.1 or a GitHub checkout) are arguments.

usage: lib_runs.py MODE OUT N_CAND N_PLAC GAP [--lib PATH] [--abs] [--first K] [--last K]
  MODE = iowa | synth
"""

import sys, json, time, datetime, logging, argparse, warnings

ap = argparse.ArgumentParser()
ap.add_argument("mode")
ap.add_argument("out")
ap.add_argument("n_cand", type=int)
ap.add_argument("n_plac", type=int)
ap.add_argument("gap", type=int)
ap.add_argument("--lib", default=None)
ap.add_argument("--abs", action="store_true")
ap.add_argument("--first", type=int, default=0)
ap.add_argument("--last", type=int, default=10**9)
ap.add_argument("--lifts", default="0.0,0.05,mde")
ap.add_argument(
    "--size_sd", type=float, default=1.0
)  # synthetic only: sd of log geo size (repo: 1.0; 0 = equal-size geos)
ap.add_argument("--pre", type=int, default=168)  # pre-period days (repo: 168)
ap.add_argument(
    "--pick", default="top"
)  # top = best-ranked design (repo); random = uniform draw among all designs the library returns
ap.add_argument("--data", default="iowa_daily_county.parquet")
a = ap.parse_args()
if a.lib:
    sys.path.insert(0, a.lib)
warnings.filterwarnings("ignore", category=FutureWarning)
import numpy as np, pandas as pd
import meridian_geox as geox
from meridian_geox import analysis as _an

LOG = []


class _H(logging.Handler):
    def emit(self, rec):
        if rec.levelno >= logging.WARNING:
            try:
                LOG.append(rec.getMessage()[:90])
            except Exception:
                LOG.append(str(rec.msg)[:90])


logging.getLogger().addHandler(_H())
logging.getLogger().setLevel(logging.WARNING)
try:
    from absl import logging as absl_logging

    absl_logging.get_absl_logger().addHandler(_H())
except Exception:
    pass

PLAC = []
_orig = _an._get_placebo_masks


def _wrapped(*args, **kw):
    m = _orig(*args, **kw)
    arr = np.asarray(m)
    PLAC.append((int(arr.shape[0]), int(np.unique(arr, axis=0).shape[0])))
    return m


_an._get_placebo_masks = _wrapped

PRE, TEST = a.pre, 28
HAS_ABS = "use_absolute_placebo_statistics" in geox.AnalysisConfig.__dataclass_fields__
if a.abs and not HAS_ABS:
    sys.exit("this library version has no use_absolute_placebo_statistics")


def panels():
    if a.mode == "iowa":
        D = pd.read_parquet(a.data)
        D["location"] = D.county.astype(int).astype(str)
        full = D.pivot_table(index="date", columns="location", values="sales", aggfunc="sum")
        full = full.reindex(pd.date_range(full.index.min(), full.index.max(), freq="D")).fillna(0.0)
        starts = list(range(0, len(full) - PRE - TEST, 60))
        for i, s in enumerate(starts):
            win = full.iloc[s : s + PRE + TEST]
            keep = win.columns[(win.iloc[:PRE] > 0).mean() > 0.2]
            win = win[keep]
            long = (
                win.stack().rename("conversions").reset_index().rename(columns={"level_0": "date"})
            )
            long["spend"] = long.conversions * 0.01
            meta = dict(
                window=i,
                start=str(win.index[0].date()),
                test_start=str(win.index[PRE].date()),
                n_geos=len(keep),
                zero_cell_share=float((win.iloc[:PRE] == 0).values.mean()),
                zero_day_share=float((win.iloc[:PRE].sum(axis=1) == 0).mean()),
            )
            yield i, long, win.index, meta
    else:
        G = 50
        for w in range(100):
            rng = np.random.default_rng(500 + w)
            T = PRE + TEST
            size = np.exp(9 + a.size_sd * rng.normal(0, 1.0, G))
            t = np.arange(T)
            common = 0.15 * np.sin(2 * np.pi * t / 7) + np.cumsum(rng.normal(0, 0.01, T))
            e = np.zeros((G, T))
            eps = rng.normal(0, 0.10, (G, T))
            for k in range(1, T):
                e[:, k] = 0.5 * e[:, k - 1] + eps[:, k]
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
            yield w, long, dates, dict(window=w, n_geos=G)


res = []
for i, long, idx, meta in panels():
    if i < a.first or i > a.last:
        continue
    t0 = time.time()
    LOG.clear()
    PLAC.clear()
    row = dict(meta, gap=a.gap, n_cand=a.n_cand, n_plac=a.n_plac, abs=a.abs)
    try:
        pre = long[long.date < idx[PRE - a.gap]]
        cfg = geox.DesignConfig(
            experiment_duration=datetime.timedelta(days=TEST),
            experiment_types=geox.ExperimentType.HOLDBACK,
            methodology=geox.Methodology.TBR,
            geo_assignment_rule=geox.GeoAssignmentRule.STRATIFIED_SAMPLING,
            cell_count=1,
            design_output_count=(1 if a.pick == "top" else 100),
            n_candidates=a.n_cand,
            seed=i,
        )
        ds = geox.run_design(pre, cfg, geox.Constraints(max_conversions_percent=0.3))
        _all = list(ds.designs.values())
        _k = 0 if a.pick == "top" else int(np.random.default_rng(90000 + i).integers(0, len(_all)))
        design = _all[_k]
        row.update(n_designs_returned=len(_all), picked_rank=_k)
        cell = list(design.designs.values())[0]
        treated = sorted(cell.treatment_geos)
        mde = float(cell.minimum_detectable_effect)
        row.update(
            n_treated=len(treated),
            mde=mde,
            design_warnings=sorted(set(LOG)),
            n_excluded_dates=len(design.excluded_dates or []),
            n_excluded_geos=len(design.excluded_geos or []),
        )
        LOG.clear()
        for name, lift in [
            (n_, {"0.0": 0.0, "0.05": 0.05, "mde": mde}[n_]) for n_ in a.lifts.split(",")
        ]:
            ana = long.copy()
            m = (ana.date >= idx[PRE]) & ana.location.isin(treated)
            ana.loc[m, "conversions"] *= 1 + lift
            kw = dict(
                design=design,
                analysis_start_date=idx[PRE],
                analysis_end_date=idx[-1],
                n_placebo_candidates=a.n_plac,
            )
            if a.abs:
                kw["use_absolute_placebo_statistics"] = True
            r = geox.analyze(ana, geox.AnalysisConfig(**kw))
            m_ = r.results["cell_1"]
            for k in ["lift", "percent_lift"]:
                e = getattr(m_, k)
                row[f"{k}_{name}"] = dict(
                    point=float(e.point_estimate),
                    lo=float(e.lower_bound),
                    hi=float(e.upper_bound),
                    sd=float(e.standard_deviation),
                    p=float(e.p_value),
                )
        row.update(
            analysis_warnings=sorted(set(LOG)),
            n_placebos=PLAC[0][0] if PLAC else None,
            n_distinct_placebos=PLAC[0][1] if PLAC else None,
        )
        row["sec"] = round(time.time() - t0, 1)
    except Exception as ex:
        row.update(error=repr(ex)[:400], design_warnings=sorted(set(LOG)))
    res.append(row)
    json.dump(res, open(a.out, "w"))
    print(
        i,
        row.get("sec"),
        row.get("n_placebos"),
        row.get("n_distinct_placebos"),
        row.get("error", "")[:80],
        flush=True,
    )

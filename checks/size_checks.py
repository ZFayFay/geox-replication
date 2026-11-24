"""Extra checks on my reimplementation, using the functions in geox.py unchanged.
The function definitions are taken from geox.py by executing its source up to (not including) the main study loop.

Checks (A/A tests only, 100 designs per study, same design-drawing rule as geox.py):
  1. p-value definition: the repo uses |t| (symmetric) two-sided p; Google's library 1.0.1 uses the equal-tailed
     2*min(left, right) p. Does the 22% depend on that choice? Also: how often does the 90% interval exclude zero?
  2. Dominant geo: Polk County (77). Re-run with Polk removed, and with the 3 largest counties removed.
  3. Synthetic negative control with equal-size geos (size sd = 0) versus the repo's size sd = 1.
usage: python checks/size_checks.py checks/size_checks.json [comma-separated variants]
"""

import sys, os, json, numpy as np, pandas as pd

OUT_PATH = os.path.abspath(sys.argv[1])
VARIANTS = (
    sys.argv[2].split(",")
    if len(sys.argv) > 2
    else [
        "iowa_all",
        "iowa_drop_polk",
        "iowa_drop_top3",
        "synthetic_size_sd_1",
        "synthetic_equal_size",
    ]
)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(REPO)  # geox.py reads iowa_weekly_county.parquet from the repo root
src = open(os.path.join(REPO, "geox.py")).read()
src = src[: src.index("studies = [(s, st)")]
import matplotlib

matplotlib.use("Agg")
ns = {}
exec(compile(src, "geox_defs", "exec"), ns)
W = ns["W"]
WIN = ns["WIN"]
I_TR = ns["I_TR"]
N_CAND = ns["N_CAND"]
strata, gen_masks, gmean, r2_val, placebos, tbr_stats, tbr_placebo_stats, analysis_windows = (
    ns[k]
    for k in [
        "strata",
        "gen_masks",
        "gmean",
        "r2_val",
        "placebos",
        "tbr_stats",
        "tbr_placebo_stats",
        "analysis_windows",
    ]
)


def aa_rows(Y, s, rng, n_aa):
    share = Y[I_TR].sum(0) / Y[I_TR].sum()
    labels = strata(Y[I_TR], seed=s)
    cand = gen_masks("stratified", share, labels, N_CAND, rng)
    r2 = r2_val(gmean(Y, cand), gmean(Y, 1 - cand))
    elig = np.where(r2 >= 0.8)[0]
    pick = rng.choice(elig, n_aa, replace=False) if len(elig) >= n_aa else np.argsort(-r2)[:n_aa]
    w = analysis_windows()
    out = []
    for j in pick:
        tm = cand[j]
        e, r_, _ = tbr_stats(
            gmean(Y, tm[None]), gmean(Y, 1 - tm[None]), w["fit"], w["rmse"], w["eff"]
        )
        t = e[0] / r_[0]
        row = dict(study=s)
        for name, da in [("da", True), ("rand", False)]:
            PT, PC = placebos(Y, tm, "stratified", labels, rng, da)
            tp = tbr_placebo_stats(Y, PT, PC, w)
            n = len(tp)
            row[f"p_abs_{name}"] = (1 + np.sum(np.abs(tp) >= abs(t))) / (1 + n)  # repo definition
            row[f"p_eq_{name}"] = (1 + min(2 * min(np.sum(tp <= t), np.sum(tp >= t)), n)) / (
                1 + n
            )  # library 1.0.1 definition
            qlo, qhi = np.quantile(tp, [0.05, 0.95])
            row[f"ci_excl0_{name}"] = bool(not (qlo <= t <= qhi))  # 90% interval excludes zero
            row[f"pt_outside_ci_{name}"] = bool(
                not (qlo <= 0 <= qhi)
            )  # point estimate outside its own interval
            row[f"sd_{name}"] = float(tp.std(ddof=1))
            row[f"med_{name}"] = float(np.median(tp))
        row["t"] = float(t)
        out.append(row)
    return out


def summarize(rows, key="study"):
    D = pd.DataFrame(rows)
    rb = np.random.default_rng(0)
    res = {"n_tests": len(D), "n_clusters": int(D[key].nunique())}
    for c in ["p_abs_da", "p_eq_da", "p_abs_rand", "p_eq_rand"]:
        g = (D[c] < 0.10).groupby(D[key]).mean().values
        bs = [g[rb.integers(0, len(g), len(g))].mean() for _ in range(4000)]
        res[c] = dict(
            fpr=float((D[c] < 0.10).mean()),
            ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            median_per_cluster=float(np.median(g)),
        )
    for c in ["ci_excl0_da", "ci_excl0_rand", "pt_outside_ci_da", "pt_outside_ci_rand"]:
        res[c] = float(D[c].mean())
    res["sd_t_within_cluster"] = float(D.groupby(key).t.std().mean())
    res["sd_placebo_da"] = float(D.sd_da.mean())
    res["sd_placebo_rand"] = float(D.sd_rand.mean())
    return res


OUT = {}
size_rank = W.mean().sort_values(ascending=False)
OUT["county_size"] = dict(
    top5_share={str(k): float(v) for k, v in (size_rank / size_rank.sum()).head(5).items()},
    sd_log_mean_weekly_sales=float(np.log(W.mean()[W.mean() > 0]).std()),
)
print(OUT["county_size"], flush=True)
studies = [(s, st) for s, st in enumerate(range(0, len(W) - WIN + 1, 8))]
for label, drop in [
    ("iowa_all", []),
    ("iowa_drop_polk", [77.0]),
    ("iowa_drop_top3", list(size_rank.index[:3])),
]:
    if label not in VARIANTS:
        continue
    rows = []
    for s, start in studies:
        rng = np.random.default_rng(1000 + s)
        Wd = W.iloc[start : start + WIN].drop(columns=drop)
        Wd = Wd.loc[:, ((Wd > 0) & Wd.notna()).all()]
        rows += aa_rows(Wd.to_numpy(float), s, rng, 100)
    OUT[label] = summarize(rows)
    OUT[label]["dropped"] = [float(d) for d in drop]
    print(label, json.dumps(OUT[label]), flush=True)
    json.dump(OUT, open(OUT_PATH, "w"), indent=1)
ar1 = ns["ar1"] if "ar1" in ns else None


def ar1(n, phi, sd, rng, shape=()):
    x = np.zeros(shape + (n,))
    e = rng.normal(0, sd, shape + (n,))
    for t in range(1, n):
        x[..., t] = phi * x[..., t - 1] + e[..., t]
    return x


for label, size_sd in [("synthetic_size_sd_1", 1.0), ("synthetic_equal_size", 0.0)]:
    if label not in VARIANTS:
        continue
    rows = []
    for sim in range(30):
        rng = np.random.default_rng(7000 + sim)
        Gs = 60
        size = np.exp(rng.normal(0, 1.0, Gs) * size_sd)  # same draws as the repo; scaled
        Y = (
            size[:, None]
            * np.exp(ar1(WIN, 0.8, 0.05, rng)[None, :] + ar1(WIN, 0.5, 0.08, rng, (Gs,)))
        ).T
        rows += aa_rows(Y, sim, rng, 50)
    OUT[label] = summarize(rows)
    print(label, json.dumps(OUT[label]), flush=True)
    json.dump(OUT, open(OUT_PATH, "w"), indent=1)

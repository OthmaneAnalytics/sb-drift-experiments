#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from sbdrift.estimator import DriftEstimator
from sbdrift.models import load_model_from_config
from sbdrift.utils import load_yaml

ROOT = Path(__file__).resolve().parents[1]


def summarize_run(config_path: str, model_dir: str, run_tag: str, label: str, setting: str, out_dir: Path) -> pd.DataFrame:
    ycfg = load_yaml(ROOT / config_path)
    model = load_model_from_config(ycfg)

    t0 = float(ycfg["rate_time"])
    xi0 = np.asarray(ycfg["xi0"], dtype=float)
    x0 = np.asarray(ycfg["clt_point"]["x0"], dtype=float)

    raw_dir = ROOT / "results" / "raw" / "rate" / model_dir / run_tag
    sel = pd.read_csv(raw_dir / "selected_metrics.csv")
    sel = sel[sel["method"] == "lepski"].copy().reset_index(drop=True)

    rows = []
    for i, (_, r) in enumerate(sel.iterrows(), start=1):
        M = int(r["M"])
        seed = int(r["seed"])
        rep = int(r["rep"])
        h = float(r["selected_h"])

        rng = np.random.default_rng(seed)
        xs, xu = model.sample(M, rng)
        est = DriftEstimator(model=model, xs=xs, xu=xu)
        details = est.point_details(t=t0, x=x0, xi=xi0, h=h)

        Fvals = np.asarray(details["F_matrix"], dtype=float).reshape(-1)
        w = np.asarray(details["kernel_weights"], dtype=float).reshape(-1)
        xu_arr = np.asarray(xu, dtype=float).reshape(-1)

        raw_D = Fvals * w
        raw_N = xu_arr * raw_D

        rows.append({
            "label": label,
            "setting": setting,
            "M": M,
            "rep": rep,
            "seed": seed,
            "selected_h": h,
            "sup_err": float(r["sup_err"]),
            "ise": float(r["ise"]),
            "boundary": float(r["boundary"]),
            "gap_to_oracle": float(r["gap_to_oracle"]),
            "D_hat": float(np.asarray(details["D_hat"]).reshape(-1)[0]),
            "var_raw_num": float(np.var(raw_N, ddof=1)),
            "var_raw_den": float(np.var(raw_D, ddof=1)),
        })

        if i % 20 == 0 or i == len(sel):
            print(f"[{label} {setting}] {i}/{len(sel)} done")

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / f"{label}_{setting}_replicates.csv", index=False)
    return df


def main():
    tag = Path("/tmp/sb_stress_tag.txt").read_text().strip()
    out_dir = ROOT / "results" / "processed" / "stress" / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("configs/gg_1d.yaml",                    "gg_1d",             f"gg1_stress_compact_{tag}", "GG1", "compact"),
        ("configs/gg_1d_stress_wide.yaml",        "gg_1d_wide",        f"gg1_stress_wide_{tag}",    "GG1", "wide"),
        ("configs/mm_1d.yaml",                    "mm_1d",             f"mm1_stress_compact_{tag}", "MM1", "compact"),
        ("configs/mm_1d_stress_wide_strong.yaml", "mm_1d_wide_strong", f"mm1_stress_wide_{tag}",    "MM1", "wide"),
    ]

    for config_path, model_dir, run_tag, label, setting in cases:
        summarize_run(config_path, model_dir, run_tag, label, setting, out_dir)

    summary_rows = []
    for label in ["GG1", "MM1"]:
        comp = pd.read_csv(out_dir / f"{label}_compact_replicates.csv")
        wide = pd.read_csv(out_dir / f"{label}_wide_replicates.csv")
        tau = 0.25 * comp["D_hat"].median()

        for setting, df in [("compact", comp), ("wide", wide)]:
            summary_rows.append({
                "label": label,
                "setting": setting,
                "M": int(df["M"].iloc[0]),
                "R": int(len(df)),
                "tau_rho": float(tau),
                "mean_var_raw_num": float(df["var_raw_num"].mean()),
                "mean_var_raw_den": float(df["var_raw_den"].mean()),
                "p_small_denom": float((df["D_hat"] <= tau).mean()),
                "median_sup_err": float(df["sup_err"].median()),
                "median_ise": float(df["ise"].median()),
                "mean_gap_to_oracle": float(df["gap_to_oracle"].mean()),
                "mean_selected_h": float(df["selected_h"].mean()),
                "boundary_rate": float(df["boundary"].mean()),
            })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "stress_summary.csv", index=False)

    latest = ROOT / "results" / "processed" / "stress" / "latest"
    latest.mkdir(parents=True, exist_ok=True)
    summary.to_csv(latest / "stress_summary.csv", index=False)

    def get(label, setting, col):
        return float(summary[(summary.label == label) & (summary.setting == setting)][col].iloc[0])

    values = {
        "tag": tag,
        "GG1_INFL_VARN": get("GG1", "wide", "mean_var_raw_num") / get("GG1", "compact", "mean_var_raw_num"),
        "GG1_INFL_VARD": get("GG1", "wide", "mean_var_raw_den") / get("GG1", "compact", "mean_var_raw_den"),
        "GG1_INFL_MEDSUP": get("GG1", "wide", "median_sup_err") / get("GG1", "compact", "median_sup_err"),
        "GG1_INFL_MEDISE": get("GG1", "wide", "median_ise") / get("GG1", "compact", "median_ise"),
        "MM1_INFL_VARN": get("MM1", "wide", "mean_var_raw_num") / get("MM1", "compact", "mean_var_raw_num"),
        "MM1_INFL_VARD": get("MM1", "wide", "mean_var_raw_den") / get("MM1", "compact", "mean_var_raw_den"),
        "MM1_INFL_MEDSUP": get("MM1", "wide", "median_sup_err") / get("MM1", "compact", "median_sup_err"),
        "MM1_INFL_MEDISE": get("MM1", "wide", "median_ise") / get("MM1", "compact", "median_ise"),
    }

    (out_dir / "stress_values.json").write_text(json.dumps(values, indent=2), encoding="utf-8")
    (latest / "stress_values.json").write_text(json.dumps(values, indent=2), encoding="utf-8")

    print("\n=== stress summary ===")
    print(summary.to_string(index=False))
    print(f"\nWrote {out_dir / 'stress_summary.csv'}")
    print(f"Wrote {latest / 'stress_summary.csv'}")
    print(f"Wrote {out_dir / 'stress_values.json'}")
    print(f"Wrote {latest / 'stress_values.json'}")


if __name__ == "__main__":
    main()

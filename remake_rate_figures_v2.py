from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd()
FIG = ROOT / "figures"
UPLOAD = ROOT / "paper_figure_upload"
PLOT = UPLOAD / "plot_data"

FIG.mkdir(exist_ok=True)
UPLOAD.mkdir(exist_ok=True)
PLOT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 17,
    "axes.titlesize": 20,
    "axes.labelsize": 18,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 15,
    "lines.linewidth": 2.4,
    "lines.markersize": 7,
})

def pick_run(model_id: str) -> Path:
    base = ROOT / "results" / "processed" / "rate" / model_id
    dirs = [p for p in base.iterdir() if p.is_dir() and (p / "mean_sup_err.csv").exists() and (p / "mean_ise.csv").exists()]
    if not dirs:
        raise RuntimeError(f"No processed run found for {model_id}")

    preferred = [d for d in dirs if "final_rawmax_k2" in d.name.lower()]
    if preferred:
        return max(preferred, key=lambda p: p.stat().st_mtime)

    preferred = [d for d in dirs if "final" in d.name.lower()]
    if preferred:
        return max(preferred, key=lambda p: p.stat().st_mtime)

    return max(dirs, key=lambda p: p.stat().st_mtime)

def load_metric(model_id: str, metric: str):
    run = pick_run(model_id)
    fn = "mean_sup_err.csv" if metric == "sup" else "mean_ise.csv"
    df = pd.read_csv(run / fn).copy()
    ycol = "mean_sup_err" if "mean_sup_err" in df.columns else "mean_ise"
    out = df[["M", "method", ycol]].copy()
    out.columns = ["M", "method", "y"]
    out["method"] = out["method"].str.lower()
    out = out[out["method"].isin(["oracle", "lepski"])].sort_values(["method", "M"])
    out["source_run"] = str(run)
    return out

def fmt_k(x):
    x = int(x)
    return f"{x//1000}k" if x % 1000 == 0 else str(x)

def make_fig(family, ids, metric, outname, ylabel):
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.4))
    all_tables = []

    for ax, model_id in zip(axes, ids):
        df = load_metric(model_id, metric)
        all_tables.append(df.assign(model_id=model_id))

        for method, label, marker in [
            ("oracle", "Oracle $h^\\star$", "o"),
            ("lepski", "Adaptive $\\hat h$", "s"),
        ]:
            sub = df[df["method"] == method]
            ax.plot(sub["M"], sub["y"], marker=marker, label=label)

        xs = sorted(df["M"].unique())
        ax.set_xscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels([fmt_k(x) for x in xs], rotation=0)
        ax.tick_params(axis="x", pad=8)

        dim = "1D" if model_id.endswith("_1d") else "2D"
        ax.set_title(f"{family} ({dim})", pad=10)
        ax.set_xlabel("Sample size $M$")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", frameon=True)

    fig.tight_layout()
    fig.savefig(FIG / outname, bbox_inches="tight")
    fig.savefig(UPLOAD / outname, bbox_inches="tight")
    plt.close(fig)

    pd.concat(all_tables, ignore_index=True).to_csv(PLOT / outname.replace(".pdf", "_source.csv"), index=False)

    print(f"Wrote {FIG / outname}")
    print(f"Wrote {UPLOAD / outname}")
    print(f"Wrote {PLOT / outname.replace('.pdf', '_source.csv')}")

make_fig("GG", ["gg_1d", "gg_2d"], "sup", "rates_gg.pdf", "Mean sup-grid error")
make_fig("MM", ["mm_1d", "mm_2d"], "sup", "rates_mm.pdf", "Mean sup-grid error")
make_fig("GG", ["gg_1d", "gg_2d"], "ise", "rates_gg_ise.pdf", "Mean ISE")
make_fig("MM", ["mm_1d", "mm_2d"], "ise", "rates_mm_ise.pdf", "Mean ISE")

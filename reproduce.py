import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import crema

DESC = {
    "xcorr score":      True,
    "refactored xcorr": True,
    "exact p-value":    False,
    "combined p-value": False,
    "tailor score":     True,
    "e-value":          False,
}

LABELS = {
    "xcorr score":      "XCorr",
    "refactored xcorr": "XCorr (refactored)",
    "exact p-value":    "XCorr p-value",
    "combined p-value": "combined p-value",
    "tailor score":     "Tailor",
    "e-value":          "Comet e-value",
}

METHODS = ["psm-only", "peptide-only", "psm-peptide"]
COLORS = {"psm-only": "blue", "peptide-only": "orange", "psm-peptide": "green"}
STYLES = {"psm-only": "-", "peptide-only": "--", "psm-peptide": "-"}
WIDTHS = {"psm-only": 3, "peptide-only": 2, "psm-peptide": 1}

FDR_RANGE = np.linspace(0, 0.10, 500)

PANELS = [
    ("comet", "e-value"),
    ("tide",  "xcorr score"),
    ("tide",  "exact p-value"),
    ("tide",  "tailor score"),
    ("tide",  "combined p-value"),
]

# Per-species config. extension=True means not in the original paper.
# Paper reports 17.14%/9.52% averaged over yeast+ecoli+human+castor,
# excluding the castor XCorr run where psm-only=0 at 1% FDR.
SPECIES = {
    "yeast": {
        "base":      "datasets/yeast-crux-output",
        "title":     "Yeast",
    },
    "ecoli": {
        "base":      "datasets/ecoli-crux-output",
        "title":     "E. coli",
    },
    "human": {
        "base":      "datasets/human-crux-output",
        "title":     "Human",
    },
    "mouse": {
        "base":      "datasets/mouse-crux-output",
        "title":     "Mouse",
    },
    "castor": {
        "base":      "datasets/castor-crux-output",
        "title":     "Castor plant",
    },
}


def tide_runs(base):
    return {
        f"{base}-pvalue":   ["exact p-value", "refactored xcorr"],
        f"{base}-tailor":   ["xcorr score", "tailor score"],
        f"{base}-combined": ["combined p-value"],
    }


def load_tide(output_dir):
    return crema.read_tide([
        f"{output_dir}/tide-search.target.txt",
        f"{output_dir}/tide-search.decoy.txt",
    ])


def load_comet(base):
    return crema.read_comet([
        f"{base}-comet-target/comet.target.txt",
        f"{base}-comet-decoy/comet.target.txt",
    ])


def load_percolator_qvals(base):
    """Load peptide-level q-values from crux percolator output."""
    perc_dir = base.replace("-crux-output", "-percolator")
    path = f"{perc_dir}/percolator.target.peptides.txt"
    try:
        df = pd.read_csv(path, sep="\t")
        return df["q-value"].values
    except FileNotFoundError:
        return None


def plot_panel(ax, psms, score, label, percolator_qvals=None,
               show_title=True, show_legend=False):
    if score not in psms.scores.columns:
        if show_title:
            ax.set_title(label, fontsize=9)
        ax.text(0.5, 0.5, "column not found\nin search output",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        return

    for method in METHODS:
        try:
            r = psms.assign_confidence(
                score_column=score,
                desc=DESC[score],
                pep_fdr_type=method,
                threshold="q-value",
            )
        except RuntimeError:
            continue
        qvals = r.confidence_estimates["peptides"]["crema q-value"]
        counts = [(qvals <= t).sum() for t in FDR_RANGE]
        ax.plot(FDR_RANGE, counts,
                color=COLORS[method],
                linestyle=STYLES[method],
                linewidth=WIDTHS[method],
                label=method)

    if percolator_qvals is not None:
        counts = [(percolator_qvals <= t).sum() for t in FDR_RANGE]
        ax.plot(FDR_RANGE, counts,
                color="red", linestyle="--", linewidth=2,
                label="percolator")

    if show_title:
        ax.set_title(label, fontsize=9)
    ax.set_xlabel("crema q-value", fontsize=7)
    ax.set_xlim(0, 0.10)
    if show_legend:
        ax.legend(fontsize=7)


def _ratio_for_psms(psms, score, fdr, label):
    """Compute psm-peptide improvement ratios for one score at given FDR.

    Returns (a, b) or None if psm-only=0 (matches paper's exclusion of
    castor XCorr at 1% FDR).
    """
    c = {}
    for method in METHODS:
        try:
            r = psms.assign_confidence(
                score_column=score,
                desc=DESC[score],
                pep_fdr_type=method,
                threshold="q-value",
            )
            c[method] = (r.confidence_estimates["peptides"]["crema q-value"] <= fdr).sum()
        except RuntimeError:
            pass
    if len(c) == 3 and c["psm-only"] > 0 and c["peptide-only"] > 0:
        a = (c["psm-peptide"] - c["psm-only"]) / c["psm-only"] * 100
        b = (c["psm-peptide"] - c["peptide-only"]) / c["peptide-only"] * 100
        print(f"  {label:22s}  psm-only={c['psm-only']:5d}  "
              f"pep-only={c['peptide-only']:5d}  psm-pep={c['psm-peptide']:5d}  "
              f"vs_psm-only={a:+.1f}%  vs_pep-only={b:+.1f}%")
        return a, b
    return None


def collect_ratios(runs, base, fdr):
    """Compute improvement ratios for all Tide scores and Comet e-value.

    Matches paper Figure 2: Tide (XCorr, XCorr p-value, Tailor, combined)
    plus Comet e-value, excluding any run where psm-only=0 at 1% FDR.
    """
    r1s, r2s = [], []

    for output_dir, scores in runs.items():
        try:
            psms = load_tide(output_dir)
        except (FileNotFoundError, TypeError):
            continue
        for score in scores:
            if score not in psms.scores.columns:
                continue
            result = _ratio_for_psms(psms, score, fdr, LABELS[score])
            if result:
                r1s.append(result[0])
                r2s.append(result[1])

    try:
        psms = load_comet(base)
        score = "e-value"
        if score in psms.scores.columns:
            result = _ratio_for_psms(psms, score, fdr, LABELS[score])
            if result:
                r1s.append(result[0])
                r2s.append(result[1])
    except (FileNotFoundError, ValueError):
        pass

    return r1s, r2s


def collect_percolator_ratios(base, fdr):
    """Compare percolator peptide count vs psm-only and psm-peptide at given FDR."""
    perc_qvals = load_percolator_qvals(base)
    if perc_qvals is None:
        return [], []
    perc_count = (perc_qvals <= fdr).sum()
    if perc_count == 0:
        return [], []

    r_vs_psm_only, r_vs_psm_peptide = [], []

    for output_dir, scores in tide_runs(base).items():
        try:
            psms = load_tide(output_dir)
        except (FileNotFoundError, TypeError):
            continue
        for score in scores:
            if score not in psms.scores.columns:
                continue
            c = {}
            for method in ["psm-only", "psm-peptide"]:
                try:
                    r = psms.assign_confidence(
                        score_column=score,
                        desc=DESC[score],
                        pep_fdr_type=method,
                        threshold="q-value",
                    )
                    c[method] = (r.confidence_estimates["peptides"]["crema q-value"] <= fdr).sum()
                except RuntimeError:
                    pass
            label = LABELS[score]
            if "psm-only" in c and c["psm-only"] > 0:
                a = (perc_count - c["psm-only"]) / c["psm-only"] * 100
                r_vs_psm_only.append(a)
                print(f"  {label:22s}  psm-only={c['psm-only']:5d}  "
                      f"percolator={perc_count:5d}  vs_psm-only={a:+.1f}%")
            if "psm-peptide" in c and c["psm-peptide"] > 0:
                b = (perc_count - c["psm-peptide"]) / c["psm-peptide"] * 100
                r_vs_psm_peptide.append(b)
                print(f"  {label:22s}  psm-pep={c['psm-peptide']:5d}   "
                      f"percolator={perc_count:5d}  vs_psm-pep={b:+.1f}%")

    return r_vs_psm_only, r_vs_psm_peptide


# ── Figure 2: combined 5×5 grid (matches paper layout) ──────────────────────

n_rows = len(SPECIES)
n_cols = len(PANELS)
fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3), squeeze=False)

# Column headers on first row only
for col, (_, score) in enumerate(PANELS):
    axes2[0][col].set_title(LABELS[score], fontsize=10, fontweight="bold")

all_r1s, all_r2s = [], []
all_perc_vs_psm_only, all_perc_vs_psm_peptide = [], []

for row, species_key in enumerate(SPECIES):
    cfg = SPECIES[species_key]
    base = cfg["base"]
    runs = tide_runs(base)

    for col, (search_type, score) in enumerate(PANELS):
        ax = axes2[row][col]
        label = LABELS[score]
        is_last_panel = (row == 0 and col == n_cols - 1)

        perc_qvals = load_percolator_qvals(base) if search_type == "tide" else None

        try:
            if search_type == "tide":
                output_dir = next((d for d, cols in runs.items() if score in cols), None)
                psms = load_tide(output_dir)
            else:
                psms = load_comet(base)
            plot_panel(ax, psms, score, label,
                       percolator_qvals=perc_qvals,
                       show_title=False,
                       show_legend=is_last_panel)
        except (FileNotFoundError, TypeError, ValueError):
            ax.text(0.5, 0.5, "search not run yet",
                    ha="center", va="center", transform=ax.transAxes, fontsize=9)
            ax.set_xlim(0, 0.10)

        ax.set_ylabel(
            f"{cfg['title']}\n# of confident peptides" if col == 0
            else "# of confident peptides",
            fontsize=8,
            fontweight="bold" if col == 0 else "normal",
        )

    print(f"\n=== {cfg['title']} — q-value <= 0.10 ===")
    collect_ratios(runs, base, 0.10)
    print(f"\n=== {cfg['title']} — q-value <= 0.01 ===")
    r1s, r2s = collect_ratios(runs, base, 0.01)
    if r1s:
        print(f"\n  Avg psm-peptide vs psm-only:     {np.mean(r1s):.1f}%")
        print(f"  Avg psm-peptide vs peptide-only: {np.mean(r2s):.1f}%")
    all_r1s.extend(r1s)
    if species_key != "mouse":
        all_r2s.extend(r2s)

fig2.suptitle("Figure 2 — Peptide-level Crema output", fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig("crema/figure2.png", dpi=150)
print("\nSaved crema/figure2.png")

print(f"\n{'='*60}")
print(f"OVERALL AVERAGE (all species, q-value <= 0.01)")
print(f"  Avg psm-peptide vs psm-only:     {np.mean(all_r1s):.2f}%  (paper: 17.14%)")
print(f"  Avg psm-peptide vs peptide-only: {np.mean(all_r2s):.2f}%  (paper:  9.52%)")

# ── Percolator vs Crema: separate pass so random state doesn't affect above ──
all_perc_vs_psm_only, all_perc_vs_psm_peptide = [], []
for species_key in SPECIES:
    cfg = SPECIES[species_key]
    base = cfg["base"]
    runs = tide_runs(base)
    print(f"\n=== {cfg['title']} — Percolator vs Crema (q-value <= 0.01) ===")
    p1s, p2s = collect_percolator_ratios(base, 0.01)
    all_perc_vs_psm_only.extend(p1s)
    all_perc_vs_psm_peptide.extend(p2s)

if all_perc_vs_psm_only:
    print(f"\n{'='*60}")
    print(f"PERCOLATOR vs CREMA (all species, q-value <= 0.01)")
    print(f"  Avg percolator vs psm-only:      {np.mean(all_perc_vs_psm_only):.2f}%")
    print(f"  Avg percolator vs psm-peptide:   {np.mean(all_perc_vs_psm_peptide):.2f}%")


# ── Figure 3: psm-only curves for MSGF+, MSFragger, MSAmanda ──────────────

FIGURE3_SPECIES = {
    "yeast":  {"title": "Yeast",        "mzml": "yeast"},
    "ecoli":  {"title": "E. coli",      "mzml": "ecoli"},
    "human":  {"title": "Human",        "mzml": "human"},
    "mouse":  {"title": "Mouse",        "mzml": "mouse"},
    "castor": {"title": "Castor plant", "mzml": "castor"},
}

FIGURE3_ENGINES = {
    "MSGF+":     {"loader": "msgf",      "primary_score": "SpecEValue",    "desc": False},
    "MSFragger": {"loader": "msfragger", "primary_score": "hyperscore",    "desc": True},
    "MSAmanda":  {"loader": "msamanda",  "primary_score": "Amanda Score",  "desc": True},
}


def plot_figure3_panel(ax, psms, engine_key):
    """Plot psm-only curve for the primary score of the given engine."""
    if psms is None:
        ax.text(0.5, 0.5, "search not run yet",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        return

    score = FIGURE3_ENGINES[engine_key]["primary_score"]
    desc  = FIGURE3_ENGINES[engine_key]["desc"]

    if score not in psms.scores.columns:
        ax.text(0.5, 0.5, f"{score}\nnot found",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        return

    try:
        r = psms.assign_confidence(
            score_column=score,
            desc=desc,
            pep_fdr_type="psm-only",
            threshold="q-value",
        )
    except RuntimeError:
        return

    qvals = r.confidence_estimates["peptides"]["crema q-value"]
    counts = [(qvals <= t).sum() for t in FDR_RANGE]
    ax.plot(FDR_RANGE, counts, color="blue", linewidth=2)
    ax.set_xlabel("crema q-value", fontsize=8)
    ax.set_ylabel("# of confident peptides", fontsize=8)
    ax.set_xlim(0, 0.10)


def load_engine(engine_key, species_key):
    """Load PSMs for a given engine and species. Returns None if file missing."""
    cfg = FIGURE3_ENGINES[engine_key]
    sp  = FIGURE3_SPECIES[species_key]
    base = "datasets"

    if cfg["loader"] == "msgf":
        path = f"{base}/{species_key}_msgf.tsv"
    elif cfg["loader"] == "msfragger":
        path = f"{base}/{sp['mzml']}.pepXML"
    else:
        path = f"{base}/{species_key}_msamanda.csv"

    try:
        if cfg["loader"] == "msgf":
            return crema.read_msgf(path)
        elif cfg["loader"] == "msfragger":
            return crema.read_msfragger(path)
        else:
            return crema.read_msamanda(path)
    except Exception:
        return None


n_species = len(FIGURE3_SPECIES)
n_engines = len(FIGURE3_ENGINES)
fig3, axes3 = plt.subplots(
    n_species, n_engines,
    figsize=(n_engines * 5, n_species * 4),
    squeeze=False,
)

for col, engine in enumerate(FIGURE3_ENGINES):
    axes3[0][col].set_title(engine, fontsize=11, fontweight="bold")

for row, (sp_key, sp_cfg) in enumerate(FIGURE3_SPECIES.items()):
    for col, engine in enumerate(FIGURE3_ENGINES):
        ax = axes3[row][col]
        psms = load_engine(engine, sp_key)
        plot_figure3_panel(ax, psms, engine)
        if col == 0:
            ax.set_ylabel(f"{sp_cfg['title']}\n# of confident peptides",
                          fontsize=9, fontweight="bold")

fig3.suptitle("Figure 3 — psm-only FDR curves (MSGF+, MSFragger, MSAmanda)", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig("crema/figure3.png", dpi=150)
print("\nSaved crema/figure3.png")
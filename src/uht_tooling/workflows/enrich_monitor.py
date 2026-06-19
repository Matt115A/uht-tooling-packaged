#!/usr/bin/env python3
"""Enrichment monitor: quantify variant enrichment from Nanopore long-read data."""

from __future__ import annotations

import csv
import logging
import math
import glob as _glob
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from Bio import SeqIO
from scipy.stats import fisher_exact

from uht_tooling.workflows.mut_rate import (
    _ensure_workspace,
    _safe_rmtree,
    run_minimap2,
    setup_logging,
)


# ── FASTQ helpers ─────────────────────────────────────────────────────────────

def expand_fastq_inputs(inputs: Iterable[str]) -> List[Path]:
    """Expand globs, deduplicate, return resolved Paths."""
    paths: List[Path] = []
    for item in inputs:
        if any(ch in item for ch in "*?[]"):
            paths.extend(Path(p) for p in _glob.glob(item))
        else:
            paths.append(Path(item))
    seen: set = set()
    unique: List[Path] = []
    for p in paths:
        r = p.resolve()
        if r not in seen:
            seen.add(r)
            unique.append(p)
    return unique


def _fastq_stem(path: Path) -> str:
    """Strip all FASTQ suffixes to get a clean sample name."""
    name = path.name
    for suffix in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


# ── Statistics ────────────────────────────────────────────────────────────────

def _wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Wilson score 95 % CI for a proportion."""
    if n == 0:
        return 0.0, 1.0
    phat = k / n
    denom = 1 + z ** 2 / n
    center = (phat + z ** 2 / (2 * n)) / denom
    half = z * math.sqrt((phat * (1 - phat) + z ** 2 / (4 * n)) / n) / denom
    return max(0.0, center - half), min(1.0, center + half)


def _enrichment_stats(
    a: int, b: int, c: int, d: int
) -> Tuple[float, float, float, float]:
    """
    Fisher exact test + Woolf (log-normal) 95 % CI for the OR.

    2×2 table:
        | a  target, sample   | b  background, sample   |
        | c  target, baseline | d  background, baseline |

    Returns: (OR, p_value, ci_lo, ci_hi)
    """
    or_val, p_val = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    # Woolf CI — 0.5 continuity correction for zero cells
    a2, b2, c2, d2 = (v + 0.5 if v == 0 else float(v) for v in (a, b, c, d))
    log_or = math.log(a2 * d2 / (b2 * c2))
    se = math.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    ci_lo = math.exp(log_or - 1.96 * se)
    ci_hi = math.exp(log_or + 1.96 * se)
    return float(or_val), float(p_val), ci_lo, ci_hi


def _bootstrap_enrichment_scores(
    n_pos_sort: int,
    n_neg_sort: int,
    n_pos_input: int,
    n_neg_input: int,
    n_bootstrap: int = 2000,
) -> Dict[str, object]:
    """
    Compute Baret et al. η and Zinchenko et al. η' with bootstrap 95 % CI.

    Baret:      η  = (N+_sort / N-_sort)   / (N+_input / N-_input)   [odds ratio]
    Zinchenko:  η' = (N+_sort / N_sort)    / (N+_input / N_input)    [relative risk]

    Both CIs are obtained by drawing n_bootstrap independent samples from
    Binomial(n_sort, p̂_sort) and Binomial(n_input, p̂_input) and taking the
    2.5th–97.5th percentiles of the resulting score distribution.
    """
    n_sort  = n_pos_sort  + n_neg_sort
    n_input = n_pos_input + n_neg_input

    # Point estimates with 0.5 continuity correction for zero cells
    def _cc(a: int, b: int) -> Tuple[float, float]:
        if a == 0 or b == 0:
            return a + 0.5, b + 0.5
        return float(a), float(b)

    ps, ns = _cc(n_pos_sort,  n_neg_sort)
    pi, ni = _cc(n_pos_input, n_neg_input)

    eta_point       = (ps / ns) / (pi / ni)
    eta_prime_point = (ps / (ps + ns)) / (pi / (pi + ni))

    # Bootstrap CIs
    rng   = np.random.default_rng(42)
    p_hat_sort  = n_pos_sort  / n_sort  if n_sort  > 0 else 0.5
    p_hat_input = n_pos_input / n_input if n_input > 0 else 0.5

    b_pos_sort  = rng.binomial(n_sort,  p_hat_sort,  n_bootstrap).astype(float)
    b_pos_input = rng.binomial(n_input, p_hat_input, n_bootstrap).astype(float)
    b_neg_sort  = n_sort  - b_pos_sort
    b_neg_input = n_input - b_pos_input

    # 0.5 correction for bootstrapped zeros
    b_pos_sort  = np.where(b_pos_sort  == 0, 0.5, b_pos_sort)
    b_neg_sort  = np.where(b_neg_sort  == 0, 0.5, b_neg_sort)
    b_pos_input = np.where(b_pos_input == 0, 0.5, b_pos_input)
    b_neg_input = np.where(b_neg_input == 0, 0.5, b_neg_input)

    eta_boot       = (b_pos_sort / b_neg_sort) / (b_pos_input / b_neg_input)
    eta_prime_boot = (b_pos_sort / (b_pos_sort + b_neg_sort)) / (b_pos_input / (b_pos_input + b_neg_input))

    return {
        "baret_eta":              round(float(eta_point),       4),
        "baret_eta_CI_95pct":     f"{np.percentile(eta_boot, 2.5):.3f}-{np.percentile(eta_boot, 97.5):.3f}",
        "zinchenko_eta_prime":    round(float(eta_prime_point), 4),
        "zinchenko_eta_prime_CI_95pct": (
            f"{np.percentile(eta_prime_boot, 2.5):.3f}-"
            f"{np.percentile(eta_prime_boot, 97.5):.3f}"
        ),
    }


# ── Alignment counting ────────────────────────────────────────────────────────

def _count_primary_alignments(sam_path: str) -> Dict[str, int]:
    """Count primary mapped reads per reference record from a SAM file."""
    counts: Dict[str, int] = {}
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            ref = read.reference_name
            if ref:
                counts[ref] = counts.get(ref, 0) + 1
    return counts


def _load_reference_ids(fasta_path: Path) -> List[str]:
    """Return ordered list of record IDs from a FASTA file."""
    return [r.id for r in SeqIO.parse(str(fasta_path), "fasta")]


# ── Figure ────────────────────────────────────────────────────────────────────

_INK   = "#1a1a1a"
_MID   = "#666666"
_LGRAY = "#cccccc"
_BLUE  = "#2166ac"
_GREEN = "#1a9850"
_CORAL = "#d73027"


def _style_ax(ax) -> None:
    ax.set_facecolor("white")
    ax.tick_params(labelsize=8.5, colors=_MID, length=3, width=0.8)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    for sp in ["bottom", "left"]:
        ax.spines[sp].set_color(_LGRAY)
        ax.spines[sp].set_linewidth(0.8)


def _superscript(n: int) -> str:
    sup = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
           "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻"}
    return "".join(sup.get(c, c) for c in str(n))


_AMBER = "#f46d43"


def _make_summary_figure(
    sample_order: List[str],
    fractions: Dict[str, float],
    ci_lo: Dict[str, float],
    ci_hi: Dict[str, float],
    or_vals: Dict[str, float],
    or_ci_lo: Dict[str, float],
    or_ci_hi: Dict[str, float],
    p_vals: Dict[str, float],
    eta_vals: Dict[str, float],
    eta_ci_lo: Dict[str, float],
    eta_ci_hi: Dict[str, float],
    eta_prime_vals: Dict[str, float],
    eta_prime_ci_lo: Dict[str, float],
    eta_prime_ci_hi: Dict[str, float],
    target_id: str,
    background_id: str,
    baseline_sample: str,
    out_path: Path,
) -> None:
    """Three-panel figure: proportional bars (A), Fisher OR forest plot (B), Baret/Zinchenko scores (C)."""
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as mticker

    n = len(sample_order)
    selected = [s for s in sample_order if s != baseline_sample]

    fig_h = max(5.5, 1.1 * n + 2.0)
    fig = plt.figure(figsize=(16, fig_h), dpi=200)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(
        1, 3, width_ratios=[2.2, 1.4, 1.4],
        left=0.13, right=0.97, bottom=0.12, top=0.90, wspace=0.45,
    )

    # ── Panel A: stacked proportion bars ──────────────────────────────────────
    axA = fig.add_subplot(gs[0])
    _style_ax(axA)

    y = np.arange(n)
    bar_h = 0.55
    baseline_frac = fractions[baseline_sample]
    fracs_ordered = [fractions[s] for s in sample_order]
    ci_lo_ord = [ci_lo[s] for s in sample_order]
    ci_hi_ord = [ci_hi[s] for s in sample_order]

    axA.barh(y, fracs_ordered, height=bar_h, color=_GREEN, label=target_id, linewidth=0)
    axA.barh(
        y,
        [1.0 - f for f in fracs_ordered],
        height=bar_h,
        left=fracs_ordered,
        color=_CORAL,
        label=background_id,
        linewidth=0,
    )
    for i, (frac, lo, hi) in enumerate(zip(fracs_ordered, ci_lo_ord, ci_hi_ord)):
        axA.plot([lo, hi], [y[i], y[i]], color="white", linewidth=2.5, zorder=4)
        axA.plot([lo, hi], [y[i], y[i]], color=_INK, linewidth=1.2, zorder=5)
        axA.plot(frac, y[i], "o", color=_INK, markersize=4.5, zorder=6)
    for i, frac in enumerate(fracs_ordered):
        axA.text(
            frac / 2, y[i], f"{frac * 100:.1f}%",
            ha="center", va="center", fontsize=8, color="white", fontweight="bold",
        )
    axA.axvline(baseline_frac, color=_INK, linestyle="--", linewidth=1.0, alpha=0.55)
    axA.text(
        baseline_frac + 0.01, n / 2,
        "unselected\nbaseline",
        ha="left", va="center", fontsize=7.5, color=_MID, style="italic",
    )
    axA.set_yticks(y)
    axA.set_yticklabels(sample_order, fontsize=9)
    axA.set_xlim(0, 1)
    axA.set_xlabel("fraction of mapped reads", fontsize=9, color=_MID)
    axA.set_title(
        f"{target_id} fraction before and after selection",
        fontsize=10.5, color=_INK, fontweight="bold", pad=8,
    )
    axA.legend(loc="lower right", fontsize=8.5, framealpha=0.9, edgecolor=_LGRAY)
    axA.text(-0.24, 1.05, "A", fontsize=16, fontweight="bold", color=_INK,
             transform=axA.transAxes)

    # ── Panel B: forest plot ───────────────────────────────────────────────────
    axB = fig.add_subplot(gs[1])
    _style_ax(axB)
    axB.spines["left"].set_visible(False)

    y2 = np.arange(len(selected))
    for i, s in enumerate(selected):
        ov, lo, hi = or_vals[s], or_ci_lo[s], or_ci_hi[s]
        p = p_vals[s]
        axB.plot([lo, hi], [y2[i], y2[i]], color=_LGRAY, linewidth=2.5, solid_capstyle="butt", zorder=2)
        axB.plot([lo, lo], [y2[i] - 0.12, y2[i] + 0.12], color=_LGRAY, linewidth=1.2, zorder=3)
        axB.plot([hi, hi], [y2[i] - 0.12, y2[i] + 0.12], color=_LGRAY, linewidth=1.2, zorder=3)
        axB.scatter([ov], [y2[i]], s=72, color=_BLUE, zorder=5, edgecolors="white", linewidth=0.8)
        if p < 1e-15:
            plabel = "p<10⁻¹⁵"
        elif p == 0.0:
            plabel = "p≈0"
        else:
            exp = math.floor(math.log10(max(p, 1e-300)))
            coef = p / 10 ** exp
            plabel = f"p={coef:.1f}×10{_superscript(exp)}"
        axB.text(hi * 1.05, y2[i] + 0.18, plabel, ha="left", va="center",
                 fontsize=7.5, color=_MID, style="italic")

    axB.axvline(1.0, color=_MID, linestyle="--", linewidth=1.0, alpha=0.65)
    axB.text(1.0, len(selected) - 0.3, "  null", ha="left", va="top",
             fontsize=7.5, color=_MID, style="italic")
    axB.set_xscale("log")

    all_los = [or_ci_lo[s] for s in selected] + [0.9]
    all_his = [or_ci_hi[s] for s in selected] + [1.1]
    xlo = min(all_los) * 0.85
    xhi = max(all_his) * 1.3
    axB.set_xlim(xlo, xhi)

    axB.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    axB.set_yticks(y2)
    axB.set_yticklabels(selected, fontsize=9)
    axB.tick_params(axis="y", length=0)
    axB.set_ylim(-0.6, len(selected) - 0.4)
    axB.set_xlabel(f"odds ratio vs. {baseline_sample} (log scale)", fontsize=9, color=_MID)
    axB.set_title(
        "Fisher exact OR\n(Woolf 95 % CI, exact p)",
        fontsize=10.0, color=_INK, fontweight="bold", pad=8,
    )
    axB.text(-0.26, 1.05, "B", fontsize=16, fontweight="bold", color=_INK,
             transform=axB.transAxes)

    # ── Panel C: Baret η and Zinchenko η' (bootstrap CI) ─────────────────────
    axC = fig.add_subplot(gs[2])
    _style_ax(axC)
    axC.spines["left"].set_visible(False)

    offset = 0.18  # vertical jitter between the two metrics per sample
    for i, s in enumerate(selected):
        # Baret η (blue, slightly above)
        ev, elo, ehi = eta_vals[s], eta_ci_lo[s], eta_ci_hi[s]
        yi = y2[i] + offset
        axC.plot([elo, ehi], [yi, yi], color=_BLUE, linewidth=2.0, solid_capstyle="butt", zorder=2, alpha=0.8)
        axC.plot([elo, elo], [yi - 0.08, yi + 0.08], color=_BLUE, linewidth=1.0, zorder=3)
        axC.plot([ehi, ehi], [yi - 0.08, yi + 0.08], color=_BLUE, linewidth=1.0, zorder=3)
        axC.scatter([ev], [yi], s=55, color=_BLUE, zorder=5, edgecolors="white", linewidth=0.7)

        # Zinchenko η' (amber, slightly below)
        epv, eplo, ephi = eta_prime_vals[s], eta_prime_ci_lo[s], eta_prime_ci_hi[s]
        yi2 = y2[i] - offset
        axC.plot([eplo, ephi], [yi2, yi2], color=_AMBER, linewidth=2.0, solid_capstyle="butt", zorder=2, alpha=0.8)
        axC.plot([eplo, eplo], [yi2 - 0.08, yi2 + 0.08], color=_AMBER, linewidth=1.0, zorder=3)
        axC.plot([ephi, ephi], [yi2 - 0.08, yi2 + 0.08], color=_AMBER, linewidth=1.0, zorder=3)
        axC.scatter([epv], [yi2], s=55, color=_AMBER, zorder=5, edgecolors="white", linewidth=0.7)

    axC.axvline(1.0, color=_MID, linestyle="--", linewidth=1.0, alpha=0.65)
    axC.set_xscale("log")

    all_c_los = [eta_ci_lo[s] for s in selected] + [eta_prime_ci_lo[s] for s in selected] + [0.9]
    all_c_his = [eta_ci_hi[s] for s in selected] + [eta_prime_ci_hi[s] for s in selected] + [1.1]
    axC.set_xlim(min(all_c_los) * 0.85, max(all_c_his) * 1.3)
    axC.get_xaxis().set_major_formatter(mticker.ScalarFormatter())

    axC.set_yticks(y2)
    axC.set_yticklabels(selected, fontsize=9)
    axC.tick_params(axis="y", length=0)
    axC.set_ylim(-0.6, len(selected) - 0.4)
    axC.set_xlabel("enrichment score (log scale)", fontsize=9, color=_MID)
    axC.set_title(
        "Baret η / Zinchenko η'\n(bootstrap 95 % CI)",
        fontsize=10.0, color=_INK, fontweight="bold", pad=8,
    )

    from matplotlib.lines import Line2D
    axC.legend(
        handles=[
            Line2D([0], [0], color=_BLUE,  marker="o", linewidth=1.5, markersize=5, label="Baret η"),
            Line2D([0], [0], color=_AMBER, marker="o", linewidth=1.5, markersize=5, label="Zinchenko η'"),
        ],
        fontsize=7.5, loc="lower right", framealpha=0.9, edgecolor=_LGRAY,
    )
    axC.text(-0.28, 1.05, "C", fontsize=16, fontweight="bold", color=_INK,
             transform=axC.transAxes)

    fig.savefig(str(out_path), facecolor="white", bbox_inches="tight")
    plt.close(fig)


# ── Main workflow ─────────────────────────────────────────────────────────────

def run_enrich_monitor(
    reference_fasta: Path,
    fastq_paths: List[Path],
    baseline_sample: str,
    output_dir: Path,
    target_id: Optional[str] = None,
    background_id: Optional[str] = None,
    work_dir: Optional[Path] = None,
    log_path: Optional[Path] = None,
    n_bootstrap: int = 2000,
) -> dict:
    """
    Quantify variant enrichment from Nanopore long-read data.

    For each FASTQ, aligns reads to the multi-sequence reference via minimap2,
    counts primary alignments per reference record, and computes enrichment
    statistics (Wilson CI on proportions, Fisher exact OR with Woolf CI) for
    each sample relative to the designated baseline (unselected pool).

    Args:
        reference_fasta: Multi-sequence FASTA file (≥2 records).
        fastq_paths: Per-sample FASTQ files (all conditions including baseline).
        baseline_sample: Sample name (FASTQ file stem) of the unselected input pool.
        output_dir: Directory where results will be written.
        target_id: FASTA record ID to treat as the enrichment target.
                   Defaults to the first record.
        background_id: FASTA record ID to treat as the background control.
                       Defaults to the second record (only if exactly 2 records).
        work_dir: Scratch directory (defaults to output_dir/tmp).
        log_path: Optional dedicated log file path.

    Returns:
        dict with keys: enrichment_stats, abundance_table, summary_figure,
        samples, baseline_sample, target_id, background_id.

    Raises:
        ValueError: if reference has <2 records, if IDs are not found, or if
                    baseline_sample does not match any FASTQ stem.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _ensure_workspace(output_dir, "write enrich-monitor outputs")

    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)

    setup_logging(str(output_dir))
    log = logging.getLogger("uht_tooling.enrich_monitor")
    if log_path:
        fh = logging.FileHandler(str(log_path))
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        log.addHandler(fh)

    # ── Scratch directory ─────────────────────────────────────────────────────
    own_work_dir = work_dir is None
    if work_dir is None:
        work_dir = output_dir / "tmp"
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # ── Reference records ─────────────────────────────────────────────────────
    ref_ids = _load_reference_ids(reference_fasta)
    if len(ref_ids) < 2:
        raise ValueError(
            f"reference_fasta must contain at least 2 records; found {len(ref_ids)}."
        )
    log.info("Reference records: %s", ref_ids)

    # ── Resolve target / background IDs ───────────────────────────────────────
    if target_id is None:
        target_id = ref_ids[0]
        log.info("Defaulting target_id to first record: %s", target_id)
    if target_id not in ref_ids:
        raise ValueError(
            f"target_id '{target_id}' not found in reference FASTA. Available: {ref_ids}"
        )

    if background_id is None:
        if len(ref_ids) == 2:
            background_id = [r for r in ref_ids if r != target_id][0]
            log.info("Defaulting background_id to second record: %s", background_id)
        else:
            raise ValueError(
                "reference_fasta has >2 records — please specify --background-id explicitly."
            )
    if background_id not in ref_ids:
        raise ValueError(
            f"background_id '{background_id}' not found in reference FASTA. Available: {ref_ids}"
        )
    if background_id == target_id:
        raise ValueError("target_id and background_id must be different records.")

    # ── Sample list ───────────────────────────────────────────────────────────
    samples: List[Dict] = [{"name": _fastq_stem(fq), "fastq": fq} for fq in fastq_paths]
    sample_names = [s["name"] for s in samples]
    if baseline_sample not in sample_names:
        raise ValueError(
            f"baseline_sample '{baseline_sample}' does not match any FASTQ stem. "
            f"Available: {sample_names}"
        )

    # ── Align and count ───────────────────────────────────────────────────────
    raw_counts: Dict[str, Dict[str, int]] = {}
    ref_fasta_str = str(Path(reference_fasta).resolve())

    for s in samples:
        name = s["name"]
        fq = str(Path(s["fastq"]).resolve())
        log.info("Aligning sample %s …", name)
        sam_path = run_minimap2(fq, ref_fasta_str, name, str(work_dir))
        counts = _count_primary_alignments(sam_path)
        raw_counts[name] = {rid: counts.get(rid, 0) for rid in ref_ids}
        log.info("  %s counts: %s", name, raw_counts[name])

    # ── Enrichment statistics ─────────────────────────────────────────────────
    baseline_target = raw_counts[baseline_sample][target_id]
    baseline_bg = raw_counts[baseline_sample][background_id]

    stats_rows: List[Dict] = []
    for s in samples:
        name = s["name"]
        tc = raw_counts[name][target_id]
        bc = raw_counts[name][background_id]
        total = tc + bc
        frac = tc / total if total > 0 else float("nan")
        wlo, whi = _wilson_ci(tc, total)

        if name == baseline_sample:
            or_val, p_val, or_lo, or_hi = 1.0, 1.0, 1.0, 1.0
            boot = {
                "baret_eta": 1.0, "baret_eta_CI_95pct": "1.000-1.000",
                "zinchenko_eta_prime": 1.0, "zinchenko_eta_prime_CI_95pct": "1.000-1.000",
            }
        else:
            or_val, p_val, or_lo, or_hi = _enrichment_stats(tc, bc, baseline_target, baseline_bg)
            boot = _bootstrap_enrichment_scores(tc, bc, baseline_target, baseline_bg, n_bootstrap)

        log2_or = math.log2(or_val) if or_val > 0 else float("-inf")

        stats_rows.append({
            "sample": name,
            target_id: tc,
            background_id: bc,
            "total_target_plus_bg": total,
            "target_fraction": round(frac, 6) if not math.isnan(frac) else "nan",
            "target_frac_CI": f"{wlo:.4f}-{whi:.4f}",
            "odds_ratio": round(or_val, 4),
            "OR_CI_95pct": f"{or_lo:.3f}-{or_hi:.3f}",
            "fisher_p": f"{p_val:.3e}",
            "log2_OR": round(log2_or, 4) if not math.isinf(log2_or) else log2_or,
            "baret_eta": boot["baret_eta"],
            "baret_eta_CI_95pct": boot["baret_eta_CI_95pct"],
            "zinchenko_eta_prime": boot["zinchenko_eta_prime"],
            "zinchenko_eta_prime_CI_95pct": boot["zinchenko_eta_prime_CI_95pct"],
            "is_baseline": name == baseline_sample,
        })

    # ── Write TSV outputs ─────────────────────────────────────────────────────
    stats_path = output_dir / "enrichment_stats.tsv"
    with stats_path.open("w", newline="") as f:
        fieldnames = list(stats_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(stats_rows)
    log.info("Wrote %s", stats_path)

    abundance_path = output_dir / "abundance_table.tsv"
    with abundance_path.open("w", newline="") as f:
        wt = csv.writer(f, delimiter="\t")
        wt.writerow(["sample"] + ref_ids)
        for s in samples:
            wt.writerow([s["name"]] + [raw_counts[s["name"]][rid] for rid in ref_ids])
    log.info("Wrote %s", abundance_path)

    sample_dirs: List[Dict] = []
    for s in samples:
        name = s["name"]
        sdir = output_dir / name
        sdir.mkdir(exist_ok=True)
        counts_path = sdir / f"{name}_counts.tsv"
        with counts_path.open("w", newline="") as f:
            wt = csv.writer(f, delimiter="\t")
            wt.writerow(["reference_id", "primary_alignments"])
            for rid in ref_ids:
                wt.writerow([rid, raw_counts[name][rid]])
        sample_dirs.append({"sample": name, "directory": str(sdir)})

    # ── Summary figure ────────────────────────────────────────────────────────
    sample_order = [s["name"] for s in samples]

    def _frac(r: Dict) -> float:
        v = r["target_fraction"]
        return float(v) if v != "nan" else 0.0

    fractions   = {r["sample"]: _frac(r) for r in stats_rows}
    ci_lo_map   = {r["sample"]: float(r["target_frac_CI"].split("-")[0]) for r in stats_rows}
    ci_hi_map   = {r["sample"]: float(r["target_frac_CI"].split("-")[1]) for r in stats_rows}
    or_map      = {r["sample"]: float(r["odds_ratio"])   for r in stats_rows}
    or_lo_map   = {r["sample"]: float(r["OR_CI_95pct"].split("-")[0]) for r in stats_rows}
    or_hi_map   = {r["sample"]: float(r["OR_CI_95pct"].split("-")[1]) for r in stats_rows}
    p_map       = {r["sample"]: float(r["fisher_p"])     for r in stats_rows}
    eta_map     = {r["sample"]: float(r["baret_eta"])    for r in stats_rows}
    eta_lo_map  = {r["sample"]: float(r["baret_eta_CI_95pct"].split("-")[0]) for r in stats_rows}
    eta_hi_map  = {r["sample"]: float(r["baret_eta_CI_95pct"].split("-")[1]) for r in stats_rows}
    etap_map    = {r["sample"]: float(r["zinchenko_eta_prime"]) for r in stats_rows}
    etap_lo_map = {r["sample"]: float(r["zinchenko_eta_prime_CI_95pct"].split("-")[0]) for r in stats_rows}
    etap_hi_map = {r["sample"]: float(r["zinchenko_eta_prime_CI_95pct"].split("-")[1]) for r in stats_rows}

    figure_path = output_dir / "enrich_summary.png"
    try:
        _make_summary_figure(
            sample_order=sample_order,
            fractions=fractions,
            ci_lo=ci_lo_map,
            ci_hi=ci_hi_map,
            or_vals=or_map,
            or_ci_lo=or_lo_map,
            or_ci_hi=or_hi_map,
            p_vals=p_map,
            eta_vals=eta_map,
            eta_ci_lo=eta_lo_map,
            eta_ci_hi=eta_hi_map,
            eta_prime_vals=etap_map,
            eta_prime_ci_lo=etap_lo_map,
            eta_prime_ci_hi=etap_hi_map,
            target_id=target_id,
            background_id=background_id,
            baseline_sample=baseline_sample,
            out_path=figure_path,
        )
        log.info("Wrote %s", figure_path)
    except Exception as exc:
        log.warning("Could not generate summary figure: %s", exc)

    # ── Cleanup scratch ───────────────────────────────────────────────────────
    if own_work_dir:
        _safe_rmtree(work_dir, allowed_base=output_dir, label="enrich_monitor tmp")

    return {
        "enrichment_stats": str(stats_path),
        "abundance_table": str(abundance_path),
        "summary_figure": str(figure_path),
        "samples": sample_dirs,
        "baseline_sample": baseline_sample,
        "target_id": target_id,
        "background_id": background_id,
    }

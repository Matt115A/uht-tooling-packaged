#!/usr/bin/env python3
"""Profiler for site-saturation mutagenesis libraries."""

from __future__ import annotations

import csv
import logging
import math
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
import pysam
from Bio.Seq import Seq

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from uht_tooling.workflows.mut_rate import (
    _ensure_workspace,
    _maybe_init_temp_workspace,
    _safe_rmtree,
    calculate_background_from_plasmid,
    compute_mismatch_stats_sam,
    describe_circular_interval,
    load_single_sequence,
    locate_gene_of_interest,
    run_minimap2,
    setup_logging,
    z_test_two_proportions,
)

plt.style.use("ggplot")

IUPAC_CODES: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")


def expand_fastq_inputs(inputs: Iterable[str]) -> List[Path]:
    paths: List[Path] = []
    for item in inputs:
        if any(ch in item for ch in "*?[]"):
            paths.extend(Path().glob(item))
        else:
            paths.append(Path(item))
    unique_paths: List[Path] = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_paths.append(path)
    return unique_paths


def create_output_directories(results_dir: Path) -> Dict[str, Path]:
    results_dir.mkdir(parents=True, exist_ok=True)
    detailed_dir = results_dir / "detailed"
    detailed_dir.mkdir(parents=True, exist_ok=True)
    return {"results_dir": results_dir, "detailed_dir": detailed_dir}


def parse_target_sites(raw_sites: Sequence[int]) -> List[int]:
    sites = sorted(set(int(site) for site in raw_sites))
    if not sites:
        raise ValueError("Provide at least one --target-site value.")
    if any(site < 1 for site in sites):
        raise ValueError("Target sites must be positive amino-acid positions (1-based).")
    return sites


def parse_scheme_map(entries: Optional[Sequence[str]]) -> Dict[int, str]:
    scheme_map: Dict[int, str] = {}
    if not entries:
        return scheme_map
    for entry in entries:
        parts = entry.split(":", 1)
        if len(parts) != 2:
            raise ValueError(
                f"Invalid scheme mapping '{entry}'. Use the form POSITION:SCHEME, for example 45:NNK."
            )
        site_text, scheme_text = parts
        site = int(site_text.strip())
        scheme = scheme_text.strip().upper()
        if site < 1:
            raise ValueError(f"Scheme site must be positive: {entry}")
        if not scheme:
            raise ValueError(f"Scheme cannot be empty: {entry}")
        invalid_chars = sorted(set(scheme) - set(IUPAC_CODES))
        if invalid_chars:
            raise ValueError(f"Unsupported IUPAC base(s) in scheme '{scheme}': {', '.join(invalid_chars)}")
        if len(scheme) != 3:
            raise ValueError(f"Scheme '{scheme}' must be exactly 3 nt long.")
        scheme_map[site] = scheme
    return scheme_map


def validate_coding_roi(hit_seq: str) -> str:
    seq_upper = hit_seq.upper()
    if len(seq_upper) % 3 != 0:
        raise ValueError("The ROI sequence must be a coding sequence with length divisible by 3.")
    protein = str(Seq(seq_upper).translate(to_stop=False))
    if "*" in protein[:-1]:
        raise ValueError("The ROI contains an internal stop codon and does not look like a clean CDS.")
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def validate_target_sites(target_sites: Sequence[int], protein_seq: str) -> None:
    protein_length = len(protein_seq)
    invalid = [site for site in target_sites if site > protein_length]
    if invalid:
        raise ValueError(
            f"Target site(s) outside the ROI protein length ({protein_length} aa): {', '.join(map(str, invalid))}"
        )


def codon_span_for_site(site: int) -> Tuple[int, int]:
    start = (site - 1) * 3
    return start, start + 3


def enumerate_degenerate_codons(scheme: str) -> List[str]:
    codons = [""]
    for char in scheme.upper():
        bases = IUPAC_CODES[char]
        codons = [prefix + base for prefix in codons for base in bases]
    return codons


def theoretical_aa_distribution(scheme: str) -> Dict[str, float]:
    codons = enumerate_degenerate_codons(scheme)
    counts = Counter(str(Seq(codon).translate()) for codon in codons)
    total = float(sum(counts.values()))
    return {aa: counts[aa] / total for aa in sorted(counts)}


def _aligned_base_map(read: pysam.AlignedSegment, ref_len: int) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    if read.is_unmapped or read.query_sequence is None:
        return mapping
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if read_pos is None or ref_pos is None:
            continue
        if not 0 <= ref_pos < ref_len:
            continue
        base = read.query_sequence[read_pos].upper()
        if base in {"A", "C", "G", "T"}:
            mapping[ref_pos] = base
    return mapping


def collect_target_site_distributions(
    sam_hit: Path,
    hit_seq: str,
    protein_seq: str,
    target_sites: Sequence[int],
) -> Dict[str, object]:
    site_results: Dict[int, Dict[str, object]] = {}
    target_positions = set()
    for site in target_sites:
        codon_start, codon_end = codon_span_for_site(site)
        target_positions.update(range(codon_start, codon_end))
        ref_codon = hit_seq[codon_start:codon_end].upper()
        site_results[site] = {
            "site": site,
            "ref_codon": ref_codon,
            "ref_aa": protein_seq[site - 1],
            "codon_counts": Counter(),
            "aa_counts": Counter(),
            "coverage": 0,
            "mutated_aa_reads": 0,
        }

    full_target_read_count = 0
    full_target_mut_loads: List[int] = []

    with pysam.AlignmentFile(str(sam_hit), "r") as samfile:
        for read in samfile.fetch():
            base_map = _aligned_base_map(read, len(hit_seq))
            per_read_mutations = 0
            covered_all_sites = True

            for site in target_sites:
                codon_start, codon_end = codon_span_for_site(site)
                codon_positions = range(codon_start, codon_end)
                if not all(pos in base_map for pos in codon_positions):
                    covered_all_sites = False
                    continue

                observed_codon = "".join(base_map[pos] for pos in codon_positions)
                observed_aa = str(Seq(observed_codon).translate())
                site_entry = site_results[site]
                site_entry["coverage"] += 1
                site_entry["codon_counts"][observed_codon] += 1
                site_entry["aa_counts"][observed_aa] += 1
                if observed_aa != site_entry["ref_aa"]:
                    site_entry["mutated_aa_reads"] += 1
                    per_read_mutations += 1

            if covered_all_sites:
                full_target_read_count += 1
                full_target_mut_loads.append(per_read_mutations)

    return {
        "sites": site_results,
        "full_target_read_count": full_target_read_count,
        "full_target_mut_loads": full_target_mut_loads,
        "target_nt_positions": sorted(target_positions),
    }


def summarize_non_target_roi(
    hit_info: Dict[str, object],
    bg_mis: int,
    bg_cov: int,
    target_positions: Sequence[int],
) -> Dict[str, object]:
    target_set = set(target_positions)
    cov = hit_info["cov"]
    mismatch = hit_info["mismatch"]
    pos_rates = hit_info["pos_rates"]

    non_target_rows: List[Dict[str, object]] = []
    total_cov = 0
    total_mis = 0
    hotspot_count = 0
    bg_rate = (bg_mis / bg_cov) if bg_cov else 0.0
    for idx, (cov_i, mis_i, rate_i) in enumerate(zip(cov, mismatch, pos_rates), start=1):
        if (idx - 1) in target_set:
            continue
        total_cov += cov_i
        total_mis += mis_i
        if cov_i == 0:
            non_target_rows.append(
                {"position": idx, "coverage": 0, "mismatches": 0, "rate": math.nan, "z_stat": 0.0, "p_value": 1.0}
            )
            continue
        z_stat, p_value = z_test_two_proportions(mis_i, cov_i, bg_mis, bg_cov)
        if rate_i > bg_rate and (p_value is None or p_value < 0.01):
            hotspot_count += 1
        non_target_rows.append(
            {
                "position": idx,
                "coverage": cov_i,
                "mismatches": mis_i,
                "rate": rate_i,
                "z_stat": z_stat,
                "p_value": p_value if p_value is not None else "",
            }
        )
    non_target_rate = (total_mis / total_cov) if total_cov else 0.0
    z_stat, p_value = z_test_two_proportions(total_mis, total_cov, bg_mis, bg_cov)
    return {
        "rows": non_target_rows,
        "coverage": total_cov,
        "mismatches": total_mis,
        "rate": non_target_rate,
        "z_stat": z_stat,
        "p_value": p_value,
        "hotspot_count": hotspot_count,
    }


def build_site_summary_rows(
    site_results: Dict[int, Dict[str, object]],
    scheme_map: Dict[int, str],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    summary_rows: List[Dict[str, object]] = []
    distribution_rows: List[Dict[str, object]] = []

    for site in sorted(site_results):
        entry = site_results[site]
        aa_counts: Counter = entry["aa_counts"]
        coverage = int(entry["coverage"])
        observed = {
            aa: aa_counts[aa] / coverage
            for aa in sorted(aa_counts)
            if coverage > 0
        }
        expected = theoretical_aa_distribution(scheme_map[site]) if site in scheme_map else {}
        all_aas = sorted(set(observed) | set(expected))
        tvd = 0.5 * sum(abs(observed.get(aa, 0.0) - expected.get(aa, 0.0)) for aa in all_aas) if expected else math.nan
        summary_rows.append(
            {
                "site": site,
                "ref_codon": entry["ref_codon"],
                "ref_aa": entry["ref_aa"],
                "scheme": scheme_map.get(site, ""),
                "codon_complete_reads": coverage,
                "mutated_aa_reads": entry["mutated_aa_reads"],
                "mutated_aa_fraction": (entry["mutated_aa_reads"] / coverage) if coverage else math.nan,
                "dominant_observed_aa": aa_counts.most_common(1)[0][0] if aa_counts else "",
                "dominant_observed_fraction": (
                    aa_counts.most_common(1)[0][1] / coverage if aa_counts and coverage else math.nan
                ),
                "observed_vs_expected_tvd": tvd,
            }
        )
        for aa in all_aas:
            distribution_rows.append(
                {
                    "site": site,
                    "ref_aa": entry["ref_aa"],
                    "scheme": scheme_map.get(site, ""),
                    "aa": aa,
                    "observed_count": aa_counts.get(aa, 0),
                    "observed_fraction": observed.get(aa, 0.0),
                    "expected_fraction": expected.get(aa, math.nan) if expected else math.nan,
                    "delta_observed_minus_expected": (
                        observed.get(aa, 0.0) - expected.get(aa, 0.0)
                        if expected
                        else math.nan
                    ),
                }
            )

    return summary_rows, distribution_rows


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _plot_summary_panels(
    results_dir: Path,
    hit_info: Dict[str, object],
    bg_rate: float,
    target_sites: Sequence[int],
    site_summary_rows: Sequence[Dict[str, object]],
    distribution_rows: Sequence[Dict[str, object]],
    non_target_rows: Sequence[Dict[str, object]],
    scheme_map: Dict[int, str],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)

    ax0 = axes[0, 0]
    positions = np.arange(1, len(hit_info["pos_rates"]) + 1)
    ax0.axhspan(0, bg_rate, color="gray", alpha=0.25, label="Plasmid background")
    ax0.plot(positions, hit_info["pos_rates"], color="#19647E", linewidth=1.5)
    for idx, site in enumerate(target_sites):
        start, _ = codon_span_for_site(site)
        ax0.axvspan(start + 1, start + 3, color="#F4D35E", alpha=0.35, label="Target codon" if idx == 0 else None)
    ax0.set_title("ROI mismatch profile")
    ax0.set_xlabel("ROI nucleotide position")
    ax0.set_ylabel("Mismatch rate")
    ax0.legend(frameon=False)

    ax1 = axes[0, 1]
    x = [row["site"] for row in site_summary_rows]
    y = [row["mutated_aa_fraction"] if not math.isnan(row["mutated_aa_fraction"]) else 0.0 for row in site_summary_rows]
    ax1.bar(x, y, color="#EE964B")
    ax1.set_title("Target-site non-reference AA fraction")
    ax1.set_xlabel("AA position")
    ax1.set_ylabel("Fraction of codon-complete reads")
    ax1.set_ylim(0, 1)

    ax2 = axes[1, 0]
    matrix = np.zeros((len(AA_ORDER), len(target_sites)), dtype=float)
    for col, site in enumerate(target_sites):
        site_rows = [row for row in distribution_rows if row["site"] == site]
        for row in site_rows:
            if row["aa"] in AA_ORDER:
                matrix[AA_ORDER.index(row["aa"]), col] = row["observed_fraction"]
    vmax = float(matrix.max()) if matrix.size else 0.0
    im = ax2.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(0.01, vmax))
    ax2.set_title("Observed AA distribution at target sites")
    ax2.set_xlabel("AA position")
    ax2.set_ylabel("Amino acid")
    ax2.set_xticks(np.arange(len(target_sites)))
    ax2.set_xticklabels(target_sites)
    ax2.set_yticks(np.arange(len(AA_ORDER)))
    ax2.set_yticklabels(AA_ORDER)
    fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)

    ax3 = axes[1, 1]
    if scheme_map:
        tvd_values = [
            row["observed_vs_expected_tvd"] if not math.isnan(row["observed_vs_expected_tvd"]) else 0.0
            for row in site_summary_rows
        ]
        ax3.bar(x, tvd_values, color="#9B5DE5")
        ax3.set_title("Observed vs expected AA distribution")
        ax3.set_xlabel("AA position")
        ax3.set_ylabel("Total variation distance")
        ax3.set_ylim(0, 1)
    else:
        xs = [row["position"] for row in non_target_rows if row["coverage"] > 0]
        ys = [row["rate"] for row in non_target_rows if row["coverage"] > 0]
        ax3.axhspan(0, bg_rate, color="gray", alpha=0.25, label="Plasmid background")
        ax3.plot(xs, ys, color="#7A306C", linewidth=1.2)
        ax3.set_title("Non-target ROI mismatch profile")
        ax3.set_xlabel("ROI nucleotide position")
        ax3.set_ylabel("Mismatch rate")
        ax3.legend(frameon=False)

    panel_png = results_dir / "summary_panels.png"
    panel_pdf = results_dir / "summary_panels.pdf"
    fig.savefig(panel_png, dpi=150)
    fig.savefig(panel_pdf)
    plt.close(fig)


def write_key_findings(
    path: Path,
    sample_name: str,
    target_sites: Sequence[int],
    full_target_read_count: int,
    full_target_mut_loads: Sequence[int],
    site_summary_rows: Sequence[Dict[str, object]],
    non_target_summary: Dict[str, object],
    bg_rate: float,
) -> None:
    mean_load = float(np.mean(full_target_mut_loads)) if full_target_mut_loads else math.nan
    wt_fraction = (
        sum(1 for value in full_target_mut_loads if value == 0) / len(full_target_mut_loads)
        if full_target_mut_loads
        else math.nan
    )
    max_site = max(
        site_summary_rows,
        key=lambda row: row["mutated_aa_fraction"] if not math.isnan(row["mutated_aa_fraction"]) else -1,
    )
    off_target_rate = non_target_summary["rate"]
    p_value = non_target_summary["p_value"]

    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"SSM profiler summary for {sample_name}\n")
        handle.write(f"{'=' * (25 + len(sample_name))}\n\n")
        handle.write(f"Target amino-acid sites: {', '.join(map(str, target_sites))}\n")
        handle.write(f"Reads spanning all target codons: {full_target_read_count}\n")
        handle.write(
            "Mean mutated target sites per fully spanning read: "
            f"{mean_load:.3f}\n" if not math.isnan(mean_load) else "Mean mutated target sites per fully spanning read: N/A\n"
        )
        handle.write(
            "Fraction fully wild-type across target sites: "
            f"{wt_fraction:.3%}\n\n" if not math.isnan(wt_fraction) else "Fraction fully wild-type across target sites: N/A\n\n"
        )
        handle.write(
            f"Most-loaded target site: {max_site['site']} ({max_site['mutated_aa_fraction']:.3%} non-reference AA)\n"
        )
        handle.write(f"Plasmid background mismatch rate: {bg_rate:.6e} per base\n")
        handle.write(
            f"Non-target ROI mismatch rate: {off_target_rate:.6e} per base "
            f"(p={p_value if p_value is not None else 'N/A'})\n"
        )
        handle.write(
            f"Non-target ROI hotspot positions above plasmid background: {non_target_summary['hotspot_count']}\n"
        )


def process_single_fastq(
    fastq_path: Path,
    region_fasta: Path,
    plasmid_fasta: Path,
    base_work_dir: Path,
    base_results_dir: Path,
    target_sites: Sequence[int],
    scheme_map: Dict[int, str],
) -> Dict[str, object]:
    _ensure_workspace(base_work_dir, "clean work directories")
    _ensure_workspace(base_results_dir, "clean result directories")

    sample_name = fastq_path.name
    if sample_name.endswith(".fastq.gz"):
        sample_name = sample_name[:-9]
    elif sample_name.endswith(".fastq"):
        sample_name = sample_name[:-6]

    work_dir = base_work_dir / sample_name
    results_dir = base_results_dir / sample_name

    if work_dir.exists():
        _safe_rmtree(work_dir, allowed_base=base_work_dir, label="work")
    if results_dir.exists():
        _safe_rmtree(results_dir, allowed_base=base_results_dir, label="results")

    work_dir.mkdir(parents=True, exist_ok=True)
    dirs = create_output_directories(results_dir)
    detailed_dir = dirs["detailed_dir"]

    setup_logging(str(results_dir))
    logging.info("--- Starting SSM analysis for sample: %s ---", sample_name)

    hit_seq, hit_id = load_single_sequence(str(region_fasta))
    plasmid_seq, _ = load_single_sequence(str(plasmid_fasta))
    protein_seq = validate_coding_roi(hit_seq)
    validate_target_sites(target_sites, protein_seq)

    roi_location = locate_gene_of_interest(plasmid_seq, hit_seq)
    logging.info(
        "ROI found at plasmid positions %s, orientation=%s",
        describe_circular_interval(roi_location["start"], roi_location["length"], len(plasmid_seq)),
        roi_location["orientation"],
    )

    sam_hit = Path(run_minimap2(str(fastq_path), str(region_fasta), "ssm_hit_alignment", str(work_dir)))
    sam_plasmid = Path(run_minimap2(str(fastq_path), str(plasmid_fasta), "ssm_plasmid_alignment", str(work_dir)))

    hit_info = compute_mismatch_stats_sam(str(sam_hit), {hit_id: hit_seq})[hit_id]
    bg_mis, bg_cov, bg_reads = calculate_background_from_plasmid(
        str(sam_plasmid),
        plasmid_seq,
        roi_location["start"],
        roi_location["length"],
    )
    bg_rate = (bg_mis / bg_cov) if bg_cov else 0.0

    target_data = collect_target_site_distributions(sam_hit, hit_seq, protein_seq, target_sites)
    site_summary_rows, distribution_rows = build_site_summary_rows(target_data["sites"], scheme_map)
    non_target_summary = summarize_non_target_roi(hit_info, bg_mis, bg_cov, target_data["target_nt_positions"])

    gene_mismatch_rows = []
    for idx, (coverage, mismatches, rate) in enumerate(
        zip(hit_info["cov"], hit_info["mismatch"], hit_info["pos_rates"]),
        start=1,
    ):
        gene_mismatch_rows.append(
            {
                "position_1based": idx,
                "coverage": coverage,
                "mismatches": mismatches,
                "mismatch_rate_per_base": "" if math.isnan(rate) else rate,
                "is_target_codon": "yes" if (idx - 1) in set(target_data["target_nt_positions"]) else "no",
            }
        )

    _write_csv(
        detailed_dir / "gene_mismatch_rates.csv",
        ["position_1based", "coverage", "mismatches", "mismatch_rate_per_base", "is_target_codon"],
        gene_mismatch_rows,
    )
    _write_csv(
        detailed_dir / "target_site_summary.csv",
        [
            "site",
            "ref_codon",
            "ref_aa",
            "scheme",
            "codon_complete_reads",
            "mutated_aa_reads",
            "mutated_aa_fraction",
            "dominant_observed_aa",
            "dominant_observed_fraction",
            "observed_vs_expected_tvd",
        ],
        site_summary_rows,
    )
    _write_csv(
        detailed_dir / "target_site_aa_distribution.csv",
        [
            "site",
            "ref_aa",
            "scheme",
            "aa",
            "observed_count",
            "observed_fraction",
            "expected_fraction",
            "delta_observed_minus_expected",
        ],
        distribution_rows,
    )
    _write_csv(
        detailed_dir / "non_target_roi_mismatch_rates.csv",
        ["position", "coverage", "mismatches", "rate", "z_stat", "p_value"],
        non_target_summary["rows"],
    )

    read_load_rows = [
        {"read_index": idx, "mutated_target_sites": value}
        for idx, value in enumerate(target_data["full_target_mut_loads"], start=1)
    ]
    _write_csv(
        detailed_dir / "target_site_load_distribution.csv",
        ["read_index", "mutated_target_sites"],
        read_load_rows,
    )

    _plot_summary_panels(
        results_dir,
        hit_info,
        bg_rate,
        target_sites,
        site_summary_rows,
        distribution_rows,
        non_target_summary["rows"],
        scheme_map,
    )

    write_key_findings(
        results_dir / "KEY_FINDINGS.txt",
        sample_name,
        target_sites,
        target_data["full_target_read_count"],
        target_data["full_target_mut_loads"],
        site_summary_rows,
        non_target_summary,
        bg_rate,
    )

    summary_path = detailed_dir / "summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Sample: {sample_name}\n")
        handle.write(f"ROI ID: {hit_id}\n")
        handle.write(f"Target sites: {', '.join(map(str, target_sites))}\n")
        handle.write(f"Reads mapped to ROI: {hit_info['mapped_reads']}\n")
        handle.write(f"Reads spanning all target sites: {target_data['full_target_read_count']}\n")
        handle.write(f"Plasmid background mismatches / coverage: {bg_mis} / {bg_cov}\n")
        handle.write(f"Plasmid background rate: {bg_rate:.6e}\n")
        handle.write(
            f"Non-target ROI mismatches / coverage: {non_target_summary['mismatches']} / {non_target_summary['coverage']}\n"
        )
        handle.write(f"Non-target ROI rate: {non_target_summary['rate']:.6e}\n")
        handle.write(
            f"Non-target ROI vs plasmid background: z={non_target_summary['z_stat']:.4f}, "
            f"p={non_target_summary['p_value'] if non_target_summary['p_value'] is not None else 'N/A'}\n"
        )

    if work_dir.exists():
        _safe_rmtree(work_dir, allowed_base=base_work_dir, label="work")

    mean_target_mut_load = (
        float(np.mean(target_data["full_target_mut_loads"]))
        if target_data["full_target_mut_loads"]
        else math.nan
    )
    return {
        "sample": sample_name,
        "results_dir": results_dir,
        "summary_path": summary_path,
        "bg_rate": bg_rate,
        "non_target_rate": non_target_summary["rate"],
        "non_target_p_value": non_target_summary["p_value"],
        "full_target_reads": target_data["full_target_read_count"],
        "mean_target_mut_load": mean_target_mut_load,
    }


def run_ssm_profiler(
    fastq_paths: Sequence[Path],
    region_fasta: Path,
    plasmid_fasta: Path,
    output_dir: Path,
    target_sites: Sequence[int],
    scheme_map: Optional[Dict[int, str]] = None,
    work_dir: Optional[Path] = None,
) -> Dict[str, object]:
    fastq_paths = [Path(path) for path in fastq_paths]
    if not fastq_paths:
        raise ValueError("No FASTQ files provided for analysis.")

    target_sites = parse_target_sites(target_sites)
    scheme_map = scheme_map or {}
    missing_scheme_sites = sorted(set(scheme_map) - set(target_sites))
    if missing_scheme_sites:
        raise ValueError(
            "Scheme(s) were provided for sites that are not in --target-site: "
            + ", ".join(map(str, missing_scheme_sites))
        )

    output_dir = Path(output_dir)
    work_dir = Path(work_dir) if work_dir is not None else output_dir / "tmp"
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    _maybe_init_temp_workspace(output_dir)
    _maybe_init_temp_workspace(work_dir)

    master_summary_path = output_dir / "master_summary.txt"
    with master_summary_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "Sample",
                    "Target_Sites",
                    "Reads_Spanning_All_Targets",
                    "Mean_Mutated_Target_Sites_Per_Read",
                    "Plasmid_Background_Rate",
                    "NonTarget_ROI_Rate",
                    "NonTarget_ROI_PValue",
                ]
            )
            + "\n"
        )

    sample_results: List[Dict[str, object]] = []
    for fastq in fastq_paths:
        result = process_single_fastq(
            fastq,
            Path(region_fasta),
            Path(plasmid_fasta),
            work_dir,
            output_dir,
            target_sites,
            scheme_map,
        )
        sample_results.append(result)
        with master_summary_path.open("a", encoding="utf-8") as handle:
            handle.write(
                "\t".join(
                    [
                        result["sample"],
                        ",".join(map(str, target_sites)),
                        str(result["full_target_reads"]),
                        "N/A" if math.isnan(result["mean_target_mut_load"]) else f"{result['mean_target_mut_load']:.6f}",
                        f"{result['bg_rate']:.6e}",
                        f"{result['non_target_rate']:.6e}",
                        str(result["non_target_p_value"] if result["non_target_p_value"] is not None else "N/A"),
                    ]
                )
                + "\n"
            )

    return {"master_summary": master_summary_path, "samples": sample_results}

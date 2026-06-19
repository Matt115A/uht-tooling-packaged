import argparse
import gzip
import logging
import os
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fuzzywuzzy import fuzz
from scipy.stats import fisher_exact, gaussian_kde
from tqdm import tqdm

DEFAULT_FLANK_MIN_RATIO = 80


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def parse_flank_spec(flanks_csv: Path) -> tuple[str, str, int, int]:
    df = pd.read_csv(flanks_csv)
    gene_start = str(df.loc[0, "gene_flanks"]).upper()
    gene_end = str(df.loc[1, "gene_flanks"]).upper()
    gene_min = int(df.loc[0, "gene_min_max"])
    gene_max = int(df.loc[1, "gene_min_max"])
    return gene_start, gene_end, gene_min, gene_max


def find_flank_matches(seq: str, flank: str, min_ratio: int) -> List[tuple[int, int, int]]:
    """Find fixed-width windows of seq whose fuzzy similarity to flank meets min_ratio.

    A literal exact-match search fails the whole read on a single base error inside the
    flank itself, which is the dominant error mode on Nanopore reads (indels, not just
    substitutions). Scoring every offset against a Levenshtein-ratio threshold tolerates
    substitutions directly, and tolerates indels immediately before/after the flank because
    the correct offset is still tried as the window slides.
    """
    flank_len = len(flank)
    if flank_len == 0 or len(seq) < flank_len:
        return []
    matches = []
    for i in range(len(seq) - flank_len + 1):
        ratio = fuzz.ratio(seq[i : i + flank_len], flank)
        if ratio >= min_ratio:
            matches.append((i, i + flank_len, ratio))
    return matches


def locate_gene_span(
    seq: str,
    gene_start: str,
    gene_end: str,
    gene_min: int,
    gene_max: int,
    min_ratio: int = DEFAULT_FLANK_MIN_RATIO,
) -> Optional[tuple[int, int, bool]]:
    """Find the (gene_span_end_of_start_flank, gene_span_start_of_end_flank, is_reverse) for seq."""
    for is_reverse, candidate in ((False, seq), (True, reverse_complement(seq))):
        start_matches = find_flank_matches(candidate, gene_start, min_ratio)
        if not start_matches:
            continue
        end_matches = find_flank_matches(candidate, gene_end, min_ratio)
        if not end_matches:
            continue
        start_matches.sort(key=lambda m: -m[2])
        end_matches.sort(key=lambda m: -m[2])
        for _, s_end, _ in start_matches[:5]:
            for e_start, _, _ in end_matches[:5]:
                if e_start <= s_end:
                    continue
                gene = candidate[s_end:e_start]
                if gene_min <= len(gene) <= gene_max:
                    return s_end, e_start, is_reverse
    return None


def extract_gene(
    seq: str,
    gene_start: str,
    gene_end: str,
    gene_min: int,
    gene_max: int,
    min_ratio: int = DEFAULT_FLANK_MIN_RATIO,
) -> Optional[str]:
    span = locate_gene_span(seq, gene_start, gene_end, gene_min, gene_max, min_ratio)
    if span is None:
        return None
    s_end, e_start, is_reverse = span
    candidate = reverse_complement(seq) if is_reverse else seq
    return candidate[s_end:e_start]


def extract_gene_and_quality(
    seq: str,
    qual: str,
    gene_start: str,
    gene_end: str,
    gene_min: int,
    gene_max: int,
    min_ratio: int = DEFAULT_FLANK_MIN_RATIO,
) -> Optional[tuple[str, str]]:
    span = locate_gene_span(seq, gene_start, gene_end, gene_min, gene_max, min_ratio)
    if span is None:
        return None
    s_end, e_start, is_reverse = span
    candidate_seq = reverse_complement(seq) if is_reverse else seq
    candidate_qual = qual[::-1] if is_reverse else qual
    return candidate_seq[s_end:e_start], candidate_qual[s_end:e_start]


def process_fastq(
    file_path: Path,
    gene_start: str,
    gene_end: str,
    gene_min: int,
    gene_max: int,
    min_ratio: int = DEFAULT_FLANK_MIN_RATIO,
) -> Dict[str, tuple[str, str]]:
    gene_reads: Dict[str, tuple[str, str]] = {}

    # Count total reads for progress bar
    total_reads = 0
    with gzip.open(file_path, "rt") as handle:
        for _ in handle:
            total_reads += 1
    total_reads = total_reads // 4

    with gzip.open(file_path, "rt") as handle:
        for _ in tqdm(range(total_reads), desc=f"Processing {file_path.name}", unit="read"):
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()  # '+'
            qual = handle.readline().strip()
            result = extract_gene_and_quality(seq, qual, gene_start, gene_end, gene_min, gene_max, min_ratio)
            if result:
                read_id = header.strip()[1:]
                gene_reads[read_id] = result
    return gene_reads


def project_quality_onto_alignment(aligned_seq: str, qual: str) -> str:
    """Re-insert gap placeholders into qual so it stays index-aligned with aligned_seq."""
    chars = []
    qi = 0
    for base in aligned_seq:
        if base == "-":
            chars.append(chr(33))  # Phred Q0 placeholder for alignment gaps
        else:
            chars.append(qual[qi])
            qi += 1
    return "".join(chars)


def align_to_reference(
    gene_reads: Dict[str, tuple[str, str]], reference: str
) -> tuple[str, Dict[str, str], Dict[str, str]]:
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_in:
        SeqIO.write(SeqRecord(Seq(reference), id="REF", description=""), tmp_in, "fasta")
        for rid, (seq, _qual) in gene_reads.items():
            SeqIO.write(SeqRecord(Seq(seq), id=rid, description=""), tmp_in, "fasta")
        tmp_in_path = tmp_in.name

    mafft = MafftCommandline(input=tmp_in_path)
    proc = subprocess.Popen(
        str(mafft),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate()

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_out:
        tmp_out.write(stdout)
        tmp_out_path = tmp_out.name

    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed with exit code {proc.returncode}:\n{stderr}")

    alignment = AlignIO.read(tmp_out_path, "fasta")

    aligned_ref = None
    aligned_reads: Dict[str, str] = {}
    for record in alignment:
        if record.id == "REF":
            aligned_ref = str(record.seq)
        else:
            aligned_reads[record.id] = str(record.seq)

    os.remove(tmp_in_path)
    os.remove(tmp_out_path)

    if aligned_ref is None:
        raise RuntimeError("Reference sequence missing from alignment output.")

    aligned_quals = {
        rid: project_quality_onto_alignment(aligned_seq, gene_reads[rid][1])
        for rid, aligned_seq in aligned_reads.items()
    }

    return aligned_ref, aligned_reads, aligned_quals


def identify_substitutions(
    ref: str,
    aligned_reads: Dict[str, str],
    aligned_quals: Optional[Dict[str, str]] = None,
    min_base_qual: int = 0,
) -> Dict[str, List[str]]:
    subs_by_read: Dict[str, List[str]] = defaultdict(list)

    aln2ref: Dict[int, Optional[int]] = {}
    ref_clean: List[str] = []
    ref_index = 0
    for aln_idx, base in enumerate(ref):
        if base != "-":
            aln2ref[aln_idx] = ref_index
            ref_clean.append(base)
            ref_index += 1
        else:
            aln2ref[aln_idx] = None
    ref_clean_seq = "".join(ref_clean)

    ref2aln: Dict[int, int] = {}
    for aln_idx, ref_idx in aln2ref.items():
        if ref_idx is not None and ref_idx not in ref2aln:
            ref2aln[ref_idx] = aln_idx

    codon_count = len(ref_clean_seq) // 3
    for read_id, seq in aligned_reads.items():
        quals = aligned_quals.get(read_id) if aligned_quals else None
        for codon_i in range(codon_count):
            start_r = codon_i * 3
            codon_ref = ref_clean_seq[start_r : start_r + 3]
            codon_read: List[str] = []
            valid = True
            diff = False
            for offset in range(3):
                ref_pos = start_r + offset
                aln_idx = ref2aln.get(ref_pos)
                if aln_idx is None:
                    valid = False
                    break
                base_q = seq[aln_idx]
                if base_q == "-":
                    valid = False
                    break
                if min_base_qual > 0 and quals is not None and (ord(quals[aln_idx]) - 33) < min_base_qual:
                    valid = False
                    break
                codon_read.append(base_q)
                if base_q != codon_ref[offset]:
                    diff = True
            if not valid or not diff:
                continue
            try:
                aa_from = str(Seq(codon_ref).translate())
                aa_to = str(Seq("".join(codon_read)).translate())
            except Exception:
                aa_from, aa_to = "?", "?"
            aa_mut = f"{aa_from}{codon_i + 1}{aa_to}"
            nt_mut = f"{codon_ref}->{''.join(codon_read)}"
            subs_by_read[read_id].append(f"{nt_mut} ({aa_mut})")

    return subs_by_read


def find_cooccurring_aa(
    subs_by_read_aa: Dict[str, List[str]],
    frequent_aa: set[str],
    output_dir: Path,
    sample_name: str,
) -> tuple[Path, Path]:
    aa_list = sorted(frequent_aa)
    aa_idx = {aa: i for i, aa in enumerate(aa_list)}

    matrix: List[List[int]] = []
    for calls in subs_by_read_aa.values():
        row = [0] * len(aa_list)
        any_selected = False
        for aa in calls:
            if aa in aa_idx:
                row[aa_idx[aa]] = 1
                any_selected = True
        if any_selected:
            matrix.append(row)

    baseline_path = output_dir / f"{sample_name}_cooccurring_AA_baseline.csv"
    fisher_path = output_dir / f"{sample_name}_cooccurring_AA_fisher.csv"

    if not matrix:
        pd.DataFrame(columns=["AA1", "AA2", "Both_Count", "AA1_Count", "AA2_Count"]).to_csv(
            baseline_path, index=False
        )
        pd.DataFrame(columns=["AA1", "AA2", "p-value"]).to_csv(fisher_path, index=False)
        return baseline_path, fisher_path

    df = pd.DataFrame(matrix, columns=aa_list)
    simple: List[tuple[str, str, int, int, int]] = []
    fisher_rows: List[tuple[str, str, float]] = []

    for i in range(len(aa_list)):
        for j in range(i + 1, len(aa_list)):
            col_a, col_b = df.iloc[:, i], df.iloc[:, j]
            both = int(((col_a == 1) & (col_b == 1)).sum())
            a_tot = int((col_a == 1).sum())
            b_tot = int((col_b == 1).sum())
            if (both >= 2) or (both > 0 and both == a_tot == b_tot):
                simple.append((aa_list[i], aa_list[j], both, a_tot, b_tot))
            table = [
                [both, int(((col_a == 1) & (col_b == 0)).sum())],
                [int(((col_a == 0) & (col_b == 1)).sum()), int(((col_a == 0) & (col_b == 0)).sum())],
            ]
            try:
                _, p_value = fisher_exact(table)
                if p_value < 0.05:
                    fisher_rows.append((aa_list[i], aa_list[j], p_value))
            except Exception:
                continue

    pd.DataFrame(simple, columns=["AA1", "AA2", "Both_Count", "AA1_Count", "AA2_Count"]).to_csv(
        baseline_path, index=False
    )
    pd.DataFrame(fisher_rows, columns=["AA1", "AA2", "p-value"]).to_csv(fisher_path, index=False)
    return baseline_path, fisher_path


def run_mutation_caller(
    template_fasta: Path,
    flanks_csv: Path,
    fastq_files: Sequence[Path],
    output_dir: Path,
    threshold: int,
    min_flank_ratio: int = DEFAULT_FLANK_MIN_RATIO,
    min_base_qual: int = 0,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> List[Dict[str, Path]]:
    template_fasta = Path(template_fasta)
    flanks_csv = Path(flanks_csv)
    fastq_files = [Path(fq) for fq in fastq_files]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.mutation_caller")
        logger.setLevel(logging.INFO)
        handler: logging.Handler
        if log_path:
            handler = logging.FileHandler(log_path, mode="w")
        else:
            handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        logger.handlers = []
        logger.addHandler(handler)
        logger.propagate = False

    try:
        if not fastq_files:
            raise ValueError("No FASTQ files provided.")

        gene_start, gene_end, gene_min, gene_max = parse_flank_spec(flanks_csv)
        template_record = next(SeqIO.parse(str(template_fasta), "fasta"))
        full_ref = str(template_record.seq)
        logger.info("Loaded template sequence of length %s.", len(full_ref))

        if full_ref.upper().startswith(gene_start) and full_ref.upper().endswith(gene_end):
            reference = full_ref[len(gene_start) : len(full_ref) - len(gene_end)]
            logger.info("Trimmed flanking regions from template.")
        else:
            reference = full_ref

        results: List[Dict[str, Path]] = []

        for fastq in tqdm(fastq_files, desc="Processing samples", unit="sample"):
            if not fastq.exists():
                logger.warning("FASTQ file %s not found; skipping.", fastq)
                continue
            sample_base = fastq.stem.replace(".fastq", "")
            sample_dir = output_dir / sample_base
            sample_dir.mkdir(parents=True, exist_ok=True)

            logger.info("Processing sample %s", sample_base)
            gene_reads = process_fastq(fastq, gene_start, gene_end, gene_min, gene_max, min_flank_ratio)
            if not gene_reads:
                logger.warning("No valid gene reads for %s; skipping.", sample_base)
                continue

            aligned_ref, aligned_reads, aligned_quals = align_to_reference(gene_reads, reference)
            substitutions = identify_substitutions(aligned_ref, aligned_reads, aligned_quals, min_base_qual)
            subs_aa = {
                rid: [item.split()[1][1:-1] for item in items if "(" in item and item.endswith(")")]
                for rid, items in substitutions.items()
            }
            counts = Counter(aa for aas in subs_aa.values() for aa in aas)
            if not counts:
                logger.warning("No amino-acid substitutions detected for %s; skipping.", sample_base)
                continue

            keys = list(counts.keys())
            values = np.array([counts[k] for k in keys], dtype=float)
            idx = np.arange(len(keys))
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            ax1.bar(idx, values)
            ax1.set_xticks(idx)
            ax1.set_xticklabels(keys, rotation=90, fontsize=8)
            ax1.set_ylabel("Count")
            ax1.set_title("Amino-Acid Substitution Frequencies")

            try:
                kde = gaussian_kde(values)
                xmin, xmax = max(1.0, values.min()), values.max()
                xs = np.logspace(np.log10(xmin), np.log10(xmax), 200)
                ax2.plot(xs, kde(xs), linewidth=2)
            except Exception:
                ax2.hist(values, bins="auto", density=True, alpha=0.6)
            ax2.set_xscale("log")
            ax2.set_xlabel("Substitution Count (log scale)")
            ax2.set_ylabel("Density")
            ax2.set_title("KDE of AA Substitution Frequencies")

            plt.tight_layout()
            plot_path = sample_dir / f"{sample_base}_aa_substitution_frequency.png"
            fig.savefig(plot_path)
            plt.close(fig)
            logger.info("Saved substitution frequency plot to %s", plot_path)

            frequent = {aa for aa, count in counts.items() if count >= threshold}
            freq_csv = sample_dir / f"{sample_base}_frequent_aa_counts.csv"
            pd.DataFrame(
                sorted(((aa, counts[aa]) for aa in frequent), key=lambda x: x[0]),
                columns=["AA", "Count"],
            ).to_csv(freq_csv, index=False)

            baseline_path, fisher_path = find_cooccurring_aa(subs_aa, frequent, sample_dir, sample_base)

            report_path = sample_dir / f"{sample_base}_report.txt"
            with report_path.open("w") as report:
                report.write(f"Sample: {sample_base}\n")
                report.write(f"Valid gene reads: {len(gene_reads)}\n")
                report.write(f"Unique AA substitutions: {len(counts)}\n")
                report.write(f"Threshold: {threshold}\n")
                report.write(f"Frequent AA substitutions (≥ {threshold}): {len(frequent)}\n\n")
                report.write("Frequent AA counts:\n")
                report.write("AA\tCount\n")
                for aa in sorted(frequent):
                    report.write(f"{aa}\t{counts[aa]}\n")
                report.write("\nGenerated files:\n")
                report.write(f"- Plot: {plot_path.name}\n")
                report.write(f"- Frequent counts: {freq_csv.name}\n")
                report.write(f"- Co-occurrence baseline: {baseline_path.name}\n")
                report.write(f"- Co-occurrence fisher:   {fisher_path.name}\n")

            results.append(
                {
                    "sample": sample_base,
                    "directory": sample_dir,
                    "plot": plot_path,
                    "frequent_counts": freq_csv,
                    "baseline": baseline_path,
                    "fisher": fisher_path,
                    "report": report_path,
                }
            )

        if not results:
            logger.warning("No outputs generated; check inputs and thresholds.")
        return results
    finally:
        if managed_logger and logger:
            for handler in list(logger.handlers):
                handler.close()
                logger.removeHandler(handler)
            logger.propagate = True


def expand_fastq_inputs(inputs: Iterable[str]) -> List[Path]:
    paths: List[Path] = []
    for item in inputs:
        if any(ch in item for ch in "*?[]"):
            paths.extend(Path().glob(item))
        else:
            paths.append(Path(item))
    unique_paths = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_paths.append(path)
    return unique_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Identify mutations from long-read sequencing without UMIs.")
    parser.add_argument("--template-fasta", required=True, type=Path, help="FASTA file containing the gene template.")
    parser.add_argument("--flanks-csv", required=True, type=Path, help="CSV describing gene flanks and length bounds.")
    parser.add_argument(
        "--fastq",
        required=True,
        nargs="+",
        help="One or more FASTQ(.gz) paths or glob patterns (e.g., data/*.fastq.gz).",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to place sample outputs.")
    parser.add_argument(
        "--threshold",
        type=int,
        default=10,
        help="Minimum AA substitution count to include in the frequent-substitution report (default: 10).",
    )
    parser.add_argument(
        "--min-flank-ratio",
        type=int,
        default=DEFAULT_FLANK_MIN_RATIO,
        help="Minimum fuzzy match ratio (0-100) for flank detection during read extraction (default: 80).",
    )
    parser.add_argument(
        "--min-base-qual",
        type=int,
        default=0,
        help=(
            "Minimum per-base Phred quality score required at every position of a codon "
            "for a substitution call at that codon (0 disables filtering, default: 0)."
        ),
    )
    parser.add_argument("--log-path", default=None, type=Path, help="Optional log file path.")
    return parser


def main(argv: Optional[Sequence[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    fastq_files = expand_fastq_inputs(args.fastq)
    run_mutation_caller(
        template_fasta=args.template_fasta,
        flanks_csv=args.flanks_csv,
        fastq_files=fastq_files,
        output_dir=args.output_dir,
        threshold=args.threshold,
        min_flank_ratio=args.min_flank_ratio,
        min_base_qual=args.min_base_qual,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

import argparse
import gzip
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fuzzywuzzy import fuzz
from tqdm import tqdm

sns.set_palette("husl")


def gc_percent(sequence: str) -> float:
    return gc_fraction(sequence) * 100


def load_probes(csv_path: Path, logger: logging.Logger) -> List[Dict[str, str]]:
    df = pd.read_csv(csv_path)
    if "upstream" not in df.columns or "downstream" not in df.columns:
        raise ValueError("Probe CSV must contain 'upstream' and 'downstream' columns.")
    probes = []
    for idx, row in df.iterrows():
        probes.append(
            {
                "upstream": str(row["upstream"]).upper(),
                "downstream": str(row["downstream"]).upper(),
                "name": row.get("name", f"probe_{idx + 1}"),
            }
        )
    logger.info("Loaded %s probe pairs from %s", len(probes), csv_path)
    return probes


def find_probe_positions(sequence: str, probe: str, min_ratio: int) -> List[tuple[int, int, int, str]]:
    positions: List[tuple[int, int, int, str]] = []
    probe_len = len(probe)
    rc_probe = str(Seq(probe).reverse_complement())
    for i in range(len(sequence) - probe_len + 1):
        window = sequence[i : i + probe_len]
        ratio_forward = fuzz.ratio(window, probe)
        if ratio_forward >= min_ratio:
            positions.append((i, i + probe_len, ratio_forward, "forward"))
        ratio_reverse = fuzz.ratio(window, rc_probe)
        if ratio_reverse >= min_ratio:
            positions.append((i, i + probe_len, ratio_reverse, "reverse"))
    return positions


def extract_inserts(
    fastq_path: Path,
    probes: List[Dict[str, str]],
    min_ratio: int,
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    inserts: List[Dict[str, object]] = []
    with gzip.open(fastq_path, "rt") as handle:
        total_reads = sum(1 for _ in handle) // 4
    logger.info("Processing %s reads from %s", total_reads, fastq_path)

    with gzip.open(fastq_path, "rt") as handle, tqdm(total=total_reads, desc="Processing reads") as progress:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()

            read_id = header.strip()[1:]
            for probe in probes:
                upstream_positions = find_probe_positions(seq, probe["upstream"], min_ratio)
                downstream_positions = find_probe_positions(seq, probe["downstream"], min_ratio)
                for up_start, up_end, up_ratio, up_strand in upstream_positions:
                    for down_start, down_end, down_ratio, down_strand in downstream_positions:
                        if up_strand == "forward" and down_strand == "forward" and down_start > up_end:
                            insert_seq = seq[up_end:down_start]
                            if insert_seq:
                                inserts.append(
                                    {
                                        "sequence": insert_seq,
                                        "length": len(insert_seq),
                                        "probe_name": probe["name"],
                                        "up_ratio": up_ratio,
                                        "down_ratio": down_ratio,
                                        "up_strand": up_strand,
                                        "down_strand": down_strand,
                                        "read_id": read_id,
                                    }
                                )
            progress.update(1)

    logger.info("Extracted %s inserts from %s", len(inserts), fastq_path)
    return inserts


def calculate_qc_metrics(inserts: List[Dict[str, object]], logger: logging.Logger) -> Dict[str, object]:
    if not inserts:
        logger.warning("No inserts found for QC analysis")
        return {}

    lengths = [insert["length"] for insert in inserts]
    sequences = [insert["sequence"] for insert in inserts]

    metrics: Dict[str, object] = {
        "total_inserts": len(inserts),
        "length_stats": {
            "mean": float(np.mean(lengths)),
            "median": float(np.median(lengths)),
            "std": float(np.std(lengths)),
            "min": int(min(lengths)),
            "max": int(max(lengths)),
            "q25": float(np.percentile(lengths, 25)),
            "q75": float(np.percentile(lengths, 75)),
        },
        "gc_content": float(np.mean([gc_percent(seq) for seq in sequences])),
        "length_distribution": Counter(lengths),
        "probe_matches": Counter(insert["probe_name"] for insert in inserts),
        "strand_combinations": Counter(f"{insert['up_strand']}-{insert['down_strand']}" for insert in inserts),
        "match_quality": {
            "up_ratio_mean": float(np.mean([insert["up_ratio"] for insert in inserts])),
            "down_ratio_mean": float(np.mean([insert["down_ratio"] for insert in inserts])),
        },
    }

    all_bases = "".join(sequences)
    base_counts = Counter(all_bases)
    total_bases = len(all_bases)
    metrics["base_composition"] = {base: count / total_bases for base, count in base_counts.items()}

    seq_counts = Counter(sequences)
    metrics["unique_sequences"] = len(seq_counts)
    metrics["duplicate_rate"] = 1 - (len(seq_counts) / len(sequences))
    logger.info("QC metrics calculated for %s inserts", len(inserts))
    return metrics


def create_qc_plots(inserts: List[Dict[str, object]], metrics: Dict[str, object], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(20, 16))

    lengths = [insert["length"] for insert in inserts]
    gc_contents = [gc_percent(insert["sequence"]) for insert in inserts]
    up_ratios = [insert["up_ratio"] for insert in inserts]
    down_ratios = [insert["down_ratio"] for insert in inserts]

    plt.subplot(3, 4, 1)
    plt.hist(lengths, bins=50, alpha=0.7, edgecolor="black")
    plt.xlabel("Insert Length (bp)")
    plt.ylabel("Frequency")
    plt.title("Insert Length Distribution")
    plt.axvline(metrics["length_stats"]["mean"], color="red", linestyle="--", label="Mean")
    plt.axvline(metrics["length_stats"]["median"], color="green", linestyle="--", label="Median")
    plt.legend()

    plt.subplot(3, 4, 2)
    plt.hist(gc_contents, bins=30, alpha=0.7, edgecolor="black")
    plt.xlabel("GC Content (%)")
    plt.ylabel("Frequency")
    plt.title("GC Content Distribution")
    plt.axvline(metrics["gc_content"], color="red", linestyle="--", label="Mean")
    plt.legend()

    plt.subplot(3, 4, 3)
    plt.scatter(up_ratios, down_ratios, alpha=0.6)
    plt.xlabel("Upstream Match Ratio (%)")
    plt.ylabel("Downstream Match Ratio (%)")
    plt.title("Probe Match Quality")

    plt.subplot(3, 4, 4)
    combo_counts = metrics["strand_combinations"]
    plt.bar(combo_counts.keys(), combo_counts.values())
    plt.xlabel("Strand Combination")
    plt.ylabel("Count")
    plt.title("Strand Combination Analysis")
    plt.xticks(rotation=45)

    plt.subplot(3, 4, 5)
    base_comp = metrics["base_composition"]
    plt.bar(base_comp.keys(), base_comp.values())
    plt.xlabel("Base")
    plt.ylabel("Frequency")
    plt.title("Base Composition")

    plt.subplot(3, 4, 6)
    probe_counts = metrics["probe_matches"]
    plt.bar(probe_counts.keys(), probe_counts.values())
    plt.xlabel("Probe")
    plt.ylabel("Insert Count")
    plt.title("Probe Performance")
    plt.xticks(rotation=45)

    plt.subplot(3, 4, 7)
    plt.scatter(lengths, gc_contents, alpha=0.6)
    plt.xlabel("Insert Length (bp)")
    plt.ylabel("GC Content (%)")
    plt.title("Length vs GC Content")

    plt.subplot(3, 4, 8)
    sorted_lengths = sorted(lengths)
    cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
    plt.plot(sorted_lengths, cumulative)
    plt.xlabel("Insert Length (bp)")
    plt.ylabel("Cumulative Fraction")
    plt.title("Cumulative Length Distribution")

    plt.subplot(3, 4, 9)
    quality_scores = [(u + d) / 2 for u, d in zip(up_ratios, down_ratios)]
    plt.hist(quality_scores, bins=30, alpha=0.7, edgecolor="black")
    plt.xlabel("Average Match Quality (%)")
    plt.ylabel("Frequency")
    plt.title("Match Quality Distribution")

    plt.subplot(3, 4, 10)
    plt.boxplot(lengths)
    plt.ylabel("Insert Length (bp)")
    plt.title("Length Statistics")

    plt.subplot(3, 4, 11)
    seq_counts = Counter(insert["sequence"] for insert in inserts)
    duplicate_counts = [count for count in seq_counts.values() if count > 1]
    if duplicate_counts:
        plt.hist(duplicate_counts, bins=20, alpha=0.7, edgecolor="black")
        plt.xlabel("Duplicate Count")
        plt.ylabel("Frequency")
        plt.title("Duplicate Distribution")
    else:
        plt.text(0.5, 0.5, "No duplicates found", ha="center", va="center", transform=plt.gca().transAxes)
        plt.title("Duplicate Distribution")

    plt.subplot(3, 4, 12)
    plt.axis("off")
    summary_text = f"""
Summary Statistics:

Total Inserts: {metrics['total_inserts']}
Mean Length: {metrics['length_stats']['mean']:.1f} bp
Median Length: {metrics['length_stats']['median']:.1f} bp
GC Content: {metrics['gc_content']:.1f}%
Unique Sequences: {metrics['unique_sequences']}
Duplicate Rate: {metrics['duplicate_rate']:.2%}
Mean Match Quality: {metrics['match_quality']['up_ratio_mean']:.1f}%
"""
    plt.text(0.05, 0.95, summary_text, transform=plt.gca().transAxes, fontsize=10, va="top", family="monospace")

    plt.tight_layout()
    plot_path = output_dir / "qc_plots.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return plot_path


def save_inserts_fasta(inserts: List[Dict[str, object]], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for idx, insert in enumerate(inserts, start=1):
            header = (
                f">{insert['probe_name']}_insert_{idx}_len{insert['length']}_"
                f"up{int(insert['up_ratio'])}_down{int(insert['down_ratio'])}"
            )
            handle.write(f"{header}\n{insert['sequence']}\n")


def save_qc_report(metrics: Dict[str, object], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        handle.write("PROFILE INSERTS QC REPORT\n")
        handle.write("=" * 50 + "\n\n")
        handle.write("SUMMARY STATISTICS:\n")
        handle.write(f"Total inserts: {metrics['total_inserts']}\n")
        handle.write(f"Mean length: {metrics['length_stats']['mean']:.2f} bp\n")
        handle.write(f"Median length: {metrics['length_stats']['median']:.2f} bp\n")
        handle.write(f"Standard deviation: {metrics['length_stats']['std']:.2f} bp\n")
        handle.write(f"Min length: {metrics['length_stats']['min']} bp\n")
        handle.write(f"Max length: {metrics['length_stats']['max']} bp\n")
        handle.write(f"GC content: {metrics['gc_content']:.2f}%\n")
        handle.write(f"Unique sequences: {metrics['unique_sequences']}\n")
        handle.write(f"Duplicate rate: {metrics['duplicate_rate']:.2%}\n\n")

        handle.write("LENGTH STATISTICS:\n")
        handle.write(f"Q25: {metrics['length_stats']['q25']:.2f} bp\n")
        handle.write(f"Q75: {metrics['length_stats']['q75']:.2f} bp\n\n")

        handle.write("PROBE PERFORMANCE:\n")
        for probe, count in metrics["probe_matches"].items():
            handle.write(f"{probe}: {count} inserts\n")
        handle.write("\n")

        handle.write("STRAND COMBINATIONS:\n")
        for combo, count in metrics["strand_combinations"].items():
            handle.write(f"{combo}: {count} inserts\n")
        handle.write("\n")

        handle.write("BASE COMPOSITION:\n")
        for base, freq in metrics["base_composition"].items():
            handle.write(f"{base}: {freq:.3f}\n")
        handle.write("\n")

        handle.write("MATCH QUALITY:\n")
        handle.write(f"Mean upstream ratio: {metrics['match_quality']['up_ratio_mean']:.2f}%\n")
        handle.write(f"Mean downstream ratio: {metrics['match_quality']['down_ratio_mean']:.2f}%\n")


def run_profile_inserts(
    probes_csv: Path,
    fastq_files: Sequence[Path],
    output_dir: Path,
    min_ratio: int = 80,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> List[Dict[str, Path]]:
    probes_csv = Path(probes_csv)
    fastq_files = [Path(path) for path in fastq_files]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.profile_inserts")
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

        probes = load_probes(probes_csv, logger)
        results: List[Dict[str, Path]] = []
        for fastq in fastq_files:
            if not fastq.exists():
                logger.warning("FASTQ file %s not found; skipping.", fastq)
                continue
            sample_base = fastq.stem.replace(".fastq", "")
            sample_dir = output_dir / sample_base
            sample_dir.mkdir(parents=True, exist_ok=True)

            inserts = extract_inserts(fastq, probes, min_ratio, logger)
            if not inserts:
                logger.warning("No inserts extracted for %s; skipping.", fastq)
                continue

            metrics = calculate_qc_metrics(inserts, logger)
            if not metrics:
                logger.warning("Metrics unavailable for %s; skipping.", fastq)
                continue

            fasta_path = sample_dir / "extracted_inserts.fasta"
            save_inserts_fasta(inserts, fasta_path)

            report_path = sample_dir / "qc_report.txt"
            save_qc_report(metrics, report_path)

            plot_path = create_qc_plots(inserts, metrics, sample_dir)

            results.append(
                {
                    "sample": sample_base,
                    "directory": sample_dir,
                    "fasta": fasta_path,
                    "report": report_path,
                    "plots": plot_path,
                }
            )

        if not results:
            logger.warning("No profile inserts outputs generated.")
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
    unique_paths: List[Path] = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_paths.append(path)
    return unique_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Profile insert sequences using probe pairs and FASTQ data.")
    parser.add_argument("--probes-csv", required=True, type=Path, help="CSV containing upstream/downstream probes.")
    parser.add_argument(
        "--fastq",
        required=True,
        nargs="+",
        help="One or more FASTQ(.gz) paths or glob patterns.",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for per-sample outputs.")
    parser.add_argument(
        "--min-ratio",
        type=int,
        default=80,
        help="Minimum fuzzy match ratio (0-100) for probe detection (default: 80).",
    )
    parser.add_argument("--log-path", default=None, type=Path, help="Optional log file path.")
    return parser


def main(argv: Optional[Sequence[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    fastq_files = expand_fastq_inputs(args.fastq)
    run_profile_inserts(
        probes_csv=args.probes_csv,
        fastq_files=fastq_files,
        output_dir=args.output_dir,
        min_ratio=args.min_ratio,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

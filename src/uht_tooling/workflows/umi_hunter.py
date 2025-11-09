import argparse
import csv
import gzip
import logging
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def load_flank_config(config_csv: Path) -> dict:
    df = pd.read_csv(config_csv)
    return {
        "umi_start": df.loc[0, "umi_flanks"],
        "umi_end": df.loc[1, "umi_flanks"],
        "umi_min": int(df.loc[0, "umi_min_max"]),
        "umi_max": int(df.loc[1, "umi_min_max"]),
        "gene_start": df.loc[0, "gene_flanks"],
        "gene_end": df.loc[1, "gene_flanks"],
    }


def build_patterns(cfg: dict) -> tuple[re.Pattern, re.Pattern]:
    pattern_umi = re.compile(
        rf"{cfg['umi_start']}([ACGT]{{{cfg['umi_min']},{cfg['umi_max']}}}){cfg['umi_end']}",
        re.IGNORECASE,
    )
    pattern_gene = re.compile(rf"{cfg['gene_start']}(.*?){cfg['gene_end']}", re.IGNORECASE)
    return pattern_umi, pattern_gene


def extract_read_info(
    seq: str,
    pattern_umi: re.Pattern,
    pattern_gene: re.Pattern,
    logger: logging.Logger,
) -> tuple[Optional[str], Optional[str]]:
    umi_match = pattern_umi.search(seq)
    gene_match = pattern_gene.search(seq)
    if umi_match and gene_match:
        return umi_match.group(1), gene_match.group(1)
    rev_seq = reverse_complement(seq)
    umi_match = pattern_umi.search(rev_seq)
    gene_match = pattern_gene.search(rev_seq)
    if umi_match and gene_match:
        return umi_match.group(1), gene_match.group(1)
    logger.debug("Failed to extract UMI/gene from read")
    return None, None


def process_fastq(
    file_path: Path,
    pattern_umi: re.Pattern,
    pattern_gene: re.Pattern,
    logger: logging.Logger,
) -> tuple[int, Dict[str, List[str]]]:
    read_count = 0
    umi_info: Dict[str, List[str]] = {}
    extracted = 0
    with gzip.open(file_path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            read_count += 1

            umi, gene = extract_read_info(seq, pattern_umi, pattern_gene, logger)
            if umi and gene:
                umi_info.setdefault(umi, []).append(gene)
                extracted += 1
            if read_count % 100000 == 0:
                logger.info("Processed %s reads so far in %s", read_count, file_path.name)
    logger.info(
        "Finished reading %s: total reads=%s, extracted pairs=%s",
        file_path,
        read_count,
        extracted,
    )
    return read_count, umi_info


def levenshtein(s1: str, s2: str) -> int:
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def percent_identity(seq1: str, seq2: str) -> float:
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 1.0
    dist = levenshtein(seq1, seq2)
    return (max_len - dist) / max_len


def cluster_umis(
    umi_info: Dict[str, List[str]],
    threshold: float,
    logger: logging.Logger,
) -> List[dict]:
    logger.info("Clustering %s unique UMIs with threshold %.2f", len(umi_info), threshold)
    sorted_umis = sorted(umi_info.items(), key=lambda item: len(item[1]), reverse=True)
    clusters: List[dict] = []
    for umi, gene_list in sorted_umis:
        count = len(gene_list)
        for cluster in clusters:
            if percent_identity(umi, cluster["rep"]) >= threshold:
                cluster["total_count"] += count
                cluster["members"][umi] = count
                cluster["gene_seqs"].extend(gene_list)
                break
        else:
            clusters.append(
                {
                    "rep": umi,
                    "total_count": count,
                    "members": {umi: count},
                    "gene_seqs": list(gene_list),
                }
            )
    logger.info("Formed %s UMI clusters", len(clusters))
    return clusters


def run_mafft_alignment(reference: SeqRecord, gene_seqs: List[str]) -> AlignIO.MultipleSeqAlignment:
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_in:
        fasta_in = tmp_in.name
        reference.id = "REF_template"
        reference.description = ""
        SeqIO.write(reference, tmp_in, "fasta")
        for i, seq in enumerate(gene_seqs):
            record = SeqRecord(Seq(seq), id=f"seq{i}", description="")
            SeqIO.write(record, tmp_in, "fasta")

    mafft_cline = MafftCommandline(input=fasta_in)
    proc = subprocess.Popen(
        str(mafft_cline),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        os.remove(fasta_in)
        raise RuntimeError(f"MAFFT failed with exit code {proc.returncode}:\n{stderr}")

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_out:
        fasta_out = tmp_out.name
        tmp_out.write(stdout)

    alignment = AlignIO.read(fasta_out, "fasta")
    os.remove(fasta_in)
    os.remove(fasta_out)
    return alignment


def generate_consensus(
    reference_record: SeqRecord,
    gene_seqs: List[str],
    mutation_threshold: float,
    logger: logging.Logger,
) -> str:
    if not gene_seqs:
        return ""
    alignment = run_mafft_alignment(reference_record, gene_seqs)
    ref_record = None
    other_records: List[SeqRecord] = []
    for record in alignment:
        if record.id == "REF_template":
            ref_record = record
        else:
            other_records.append(record)
    if ref_record is None or not other_records:
        logger.warning("Reference or read sequences missing from alignment output")
        return ""

    consensus_chars: List[str] = []
    num_reads = len(other_records)
    length = alignment.get_alignment_length()
    for idx in range(length):
        ref_base = ref_record.seq[idx]
        col_bases = [record.seq[idx] for record in other_records]
        counts = Counter(col_bases)
        most_common, count = counts.most_common(1)[0]
        freq = count / num_reads
        if most_common != ref_base and freq >= mutation_threshold:
            consensus_chars.append(most_common)
        else:
            consensus_chars.append(ref_base)
    return "".join(consensus_chars).replace("-", "")


def write_umi_csv(output_file: Path, clusters: List[dict]):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Cluster Representative", "Total Count", "Members"])
        for cluster in clusters:
            members_str = "; ".join(f"{umi}:{count}" for umi, count in cluster["members"].items())
            writer.writerow([cluster["rep"], cluster["total_count"], members_str])


def write_gene_csv(
    output_file: Path,
    clusters: List[dict],
    reference_record: SeqRecord,
    mutation_threshold: float,
    logger: logging.Logger,
) -> List[SeqRecord]:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    ungapped_ref_length = len(str(reference_record.seq).replace("-", ""))
    consensus_records: List[SeqRecord] = []
    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["Cluster Representative", "Total Count", "Consensus Gene", "Length Difference", "Members"]
        )
        clusters_to_align = [cluster for cluster in clusters if cluster["total_count"] > 0]
        for idx, cluster in enumerate(tqdm(clusters_to_align, desc="Processing UMI clusters", unit="cluster")):
            consensus = generate_consensus(reference_record, cluster["gene_seqs"], mutation_threshold, logger)
            length_diff = len(consensus) - ungapped_ref_length
            members_str = "; ".join(f"{umi}:{count}" for umi, count in cluster["members"].items())
            writer.writerow([cluster["rep"], cluster["total_count"], consensus, length_diff, members_str])
            record_id = f"{cluster['rep']}_cluster{idx + 1}"
            consensus_records.append(
                SeqRecord(Seq(consensus), id=record_id, description=f"Length diff: {length_diff}")
            )
    return consensus_records


def run_umi_hunter(
    template_fasta: Path,
    config_csv: Path,
    fastq_files: Sequence[Path],
    output_dir: Path,
    umi_identity_threshold: float = 0.9,
    consensus_mutation_threshold: float = 0.7,
    min_cluster_size: int = 1,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> List[Dict[str, Path]]:
    template_fasta = Path(template_fasta)
    config_csv = Path(config_csv)
    fastq_files = [Path(fq) for fq in fastq_files]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.umi_hunter")
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

        if min_cluster_size < 1:
            raise ValueError("Minimum cluster size must be at least 1.")

        cfg = load_flank_config(config_csv)
        pattern_umi, pattern_gene = build_patterns(cfg)
        reference_record = next(SeqIO.parse(str(template_fasta), "fasta"))

        results: List[Dict[str, Path]] = []

        for fastq in fastq_files:
            if not fastq.exists():
                logger.warning("FASTQ file %s not found; skipping.", fastq)
                continue
            sample_base = fastq.stem.replace(".fastq", "")
            sample_dir = output_dir / sample_base
            sample_dir.mkdir(parents=True, exist_ok=True)

            read_count, umi_info = process_fastq(fastq, pattern_umi, pattern_gene, logger)
            if not umi_info:
                logger.warning("No UMIs extracted for %s; skipping.", fastq)
                continue

            clusters = cluster_umis(umi_info, umi_identity_threshold, logger)
            umi_csv = sample_dir / f"{sample_base}_UMI_clusters.csv"
            write_umi_csv(umi_csv, clusters)

            significant_clusters = [
                cluster for cluster in clusters if cluster["total_count"] >= min_cluster_size
            ]
            if not significant_clusters:
                logger.info(
                    "No clusters met the minimum size threshold (%s reads) for %s.",
                    min_cluster_size,
                    sample_base,
                )

            gene_csv = sample_dir / f"{sample_base}_gene_consensus.csv"
            consensus_records = write_gene_csv(
                gene_csv,
                significant_clusters,
                reference_record,
                consensus_mutation_threshold,
                logger,
            )

            fasta_out = sample_dir / f"{sample_base}_consensuses.fasta"
            SeqIO.write(consensus_records, fasta_out, "fasta")

            results.append(
                {
                    "sample": sample_base,
                    "directory": sample_dir,
                    "umi_csv": umi_csv,
                    "gene_csv": gene_csv,
                    "fasta": fasta_out,
                    "reads": read_count,
                    "clusters": len(significant_clusters),
                    "clusters_total": len(clusters),
                }
            )

        if not results:
            logger.warning("No UMI hunter outputs generated.")
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
    parser = argparse.ArgumentParser(description="Cluster UMIs and generate consensus sequences from long-read data.")
    parser.add_argument("--template-fasta", required=True, type=Path, help="Template FASTA file.")
    parser.add_argument("--config-csv", required=True, type=Path, help="CSV describing UMI and gene flanks.")
    parser.add_argument(
        "--fastq",
        required=True,
        nargs="+",
        help="One or more FASTQ(.gz) paths or glob patterns.",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for per-sample outputs.")
    parser.add_argument(
        "--umi-identity-threshold",
        type=float,
        default=0.9,
        help="UMI clustering identity threshold (default: 0.9).",
    )
    parser.add_argument(
        "--consensus-mutation-threshold",
        type=float,
        default=0.7,
        help="Consensus mutation threshold for MAFFT-derived consensus (default: 0.7).",
    )
    parser.add_argument("--log-path", default=None, type=Path, help="Optional log file path.")
    return parser


def main(argv: Optional[Sequence[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    fastq_files = expand_fastq_inputs(args.fastq)
    run_umi_hunter(
        template_fasta=args.template_fasta,
        config_csv=args.config_csv,
        fastq_files=fastq_files,
        output_dir=args.output_dir,
        umi_identity_threshold=args.umi_identity_threshold,
        consensus_mutation_threshold=args.consensus_mutation_threshold,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

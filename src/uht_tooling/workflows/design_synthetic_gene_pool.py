"""Design pooled synthetic-gene oligos with reusable lift-out primers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from Bio.Seq import Seq

from uht_tooling.workflows.design_gene_oligos import (
    DEFAULT_CONSTANT_3PRIME_DNA,
    DEFAULT_CONSTANT_5PRIME_DNA,
    DEFAULT_STOP_CODON,
    DEFAULT_TARGET_HOST,
    MIN_OVERLAP,
    MAX_OVERLAP,
    VALID_INPUT_TYPES,
    VALID_TAG_MODES,
    _clean_sequence_text,
    _constant_edge_sequences,
    _load_codon_table,
    _load_gene_oligo_config,
    _load_tag_sequence,
    _normalize_target,
    codon_optimize_protein,
    _pick_edge_overlap,
    _read_fasta_records,
    _resolve_string,
    _setup_logger,
    calc_gc_content,
)

CELL_FREE_TARGET_HOST = "e_coli"
MAX_COMMON_PRIMER_LENGTH = 99


def _normalize_pool_protein_sequence(seq: str, target_name: str, logger) -> Tuple[str, List[str], bool, bool]:
    warnings: List[str] = []
    cleaned = _clean_sequence_text(seq)
    valid = set("ACDEFGHIKLMNPQRSTVWYX*")
    invalid = set(cleaned) - valid
    if invalid:
        raise ValueError(f"Protein input contains invalid residues: {''.join(sorted(invalid))}")
    if cleaned.endswith("*"):
        cleaned = cleaned[:-1]
        warnings.append("Stripped trailing stop symbol from protein input.")
        logger.warning("%s: %s", target_name, warnings[-1])
    if "*" in cleaned:
        raise ValueError("Protein input contains an internal stop codon.")
    if not cleaned:
        raise ValueError("Protein input only contained a stop codon.")

    added_start = False
    added_stop = False
    if cleaned.startswith("M"):
        cleaned = cleaned[1:]
        warnings.append("Moved initiating methionine onto reusable 5' constant oligo.")
        logger.warning("%s: %s", target_name, warnings[-1])
    else:
        warnings.append("Added initiating methionine on reusable 5' constant oligo.")
        logger.warning("%s: %s", target_name, warnings[-1])
        added_start = True
    warnings.append("Moved terminal stop codon onto reusable 3' constant oligo.")
    logger.warning("%s: %s", target_name, warnings[-1])
    return cleaned, warnings, added_start, added_stop


def _codon_optimize_protein_with_random_sites(
    protein_seq: str,
    target_host: str,
    codon_table: Dict[str, str],
) -> str:
    pieces: List[str] = []
    segment: List[str] = []
    for aa in protein_seq:
        if aa == "X":
            if segment:
                pieces.append(codon_optimize_protein("".join(segment), target_host=target_host, codon_table=codon_table))
                segment = []
            pieces.append("NNN")
        else:
            segment.append(aa)
    if segment:
        pieces.append(codon_optimize_protein("".join(segment), target_host=target_host, codon_table=codon_table))
    return "".join(pieces)


def _build_pool_handles(
    constant_5prime: str,
    constant_3prime: str,
    stop_codon: str,
    tag_mode: str,
    tag_seq: str,
) -> Tuple[str, str, str, str, float, float, float, float]:
    desired_5, desired_3 = _constant_edge_sequences(
        tag_mode=tag_mode,
        constant_5prime=constant_5prime,
        tag_seq=tag_seq,
        stop_codon=stop_codon,
        constant_3prime=constant_3prime,
    )
    left_len, left_tm, left_gc = _pick_edge_overlap(desired_5, use_tail=True)
    right_len, right_tm, right_gc = _pick_edge_overlap(desired_3, use_tail=False)
    handle_5 = desired_5[-left_len:]
    handle_3 = desired_3[:right_len]
    return desired_5, desired_3, handle_5, handle_3, left_tm, left_gc, right_tm, right_gc


def run_design_synthetic_gene_pool(
    sequence_fasta: Path,
    output_dir: Path,
    input_type: str = "auto",
    tag_mode: str = "none",
    log_path: Optional[Path] = None,
    logger=None,
    config: Optional[Dict[str, object]] = None,
) -> Dict[str, Path]:
    if input_type not in VALID_INPUT_TYPES:
        raise ValueError(f"input_type must be one of {sorted(VALID_INPUT_TYPES)}")
    if tag_mode not in VALID_TAG_MODES:
        raise ValueError(f"tag_mode must be one of {sorted(VALID_TAG_MODES)}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_log = Path(log_path) if log_path else output_dir / "run.log"
    managed_logger, owns_logger = _setup_logger(run_log, logger)

    try:
        records = _read_fasta_records(Path(sequence_fasta))
        gene_cfg = _load_gene_oligo_config(config)
        constant_5prime = _resolve_string(gene_cfg, "constant_5prime_dna", DEFAULT_CONSTANT_5PRIME_DNA)
        constant_3prime = _resolve_string(gene_cfg, "constant_3prime_dna", DEFAULT_CONSTANT_3PRIME_DNA)
        stop_codon = _resolve_string(gene_cfg, "stop_codon", DEFAULT_STOP_CODON)
        tag_seq = _load_tag_sequence(gene_cfg, tag_mode)
        codon_table = _load_codon_table(gene_cfg)

        desired_5, desired_3, handle_5, handle_3, left_tm, left_gc, right_tm, right_gc = _build_pool_handles(
            constant_5prime=constant_5prime,
            constant_3prime=constant_3prime,
            stop_codon=stop_codon,
            tag_mode=tag_mode,
            tag_seq=tag_seq,
        )

        pool_rows: List[Dict[str, object]] = []
        ordering_rows: List[Dict[str, object]] = []
        for record in records:
            cleaned = _clean_sequence_text(record.raw_sequence)
            resolved_input_type = input_type if input_type != "auto" else ("dna" if set(cleaned) <= {"A", "C", "G", "T"} else "protein")
            if resolved_input_type == "protein":
                protein_seq, warnings, added_start, added_stop = _normalize_pool_protein_sequence(
                    record.raw_sequence,
                    target_name=record.name,
                    logger=managed_logger,
                )
                coding_seq = _codon_optimize_protein_with_random_sites(
                    protein_seq,
                    target_host=CELL_FREE_TARGET_HOST,
                    codon_table=codon_table,
                )
                if "X" in protein_seq:
                    warnings.append("Translated protein residue X to NNN in the pooled oligo output.")
                    managed_logger.info("%s: %s", record.name, warnings[-1])
                warnings.append(
                    f"Applied host-aware codon optimization for protein input using target host '{CELL_FREE_TARGET_HOST}'."
                )
                managed_logger.info("%s: %s", record.name, warnings[-1])
            else:
                coding_seq, warnings, added_start, added_stop = _normalize_target(
                    record=record,
                    resolved_input_type=resolved_input_type,
                    tag_mode=tag_mode,
                    tag_seq=tag_seq,
                    stop_codon=stop_codon,
                    codon_table=codon_table,
                    target_host=CELL_FREE_TARGET_HOST,
                    logger=managed_logger,
                )
            pool_seq = handle_5 + coding_seq + handle_3
            pool_name = f"{record.name}_POOL"
            warning_text = " | ".join(warnings)
            pool_rows.append(
                {
                    "target_name": record.name,
                    "input_type": resolved_input_type,
                    "pool_oligo_name": pool_name,
                    "sequence_5to3": pool_seq,
                    "pool_oligo_length": len(pool_seq),
                    "added_start_codon": "yes" if added_start else "no",
                    "added_stop_codon": "yes" if added_stop else "no",
                    "warnings": warning_text,
                }
            )
            ordering_rows.append(
                {
                    "name": pool_name,
                    "sequence": pool_seq,
                    "category": "pool_oligo",
                    "target_name": record.name,
                    "order_scope": "order_per_target",
                }
            )

        forward_primer = desired_5
        reverse_primer = str(Seq(desired_3).reverse_complement())
        if len(forward_primer) > MAX_COMMON_PRIMER_LENGTH or len(reverse_primer) > MAX_COMMON_PRIMER_LENGTH:
            raise ValueError(
                "The common primer pair exceeds the 99 bp ordering limit. Shorten the constant/tag payload for this cell-free design."
            )
        common_primers = [
            {
                "primer_name": "POOL_CONST_F",
                "sequence_5to3": forward_primer,
                "orientation": "forward",
                "adds_sequence": desired_5,
                "anneal_sequence": handle_5,
                "anneal_length_nt": len(handle_5),
                "anneal_tm_c": f"{left_tm:.2f}",
                "anneal_gc_percent": f"{left_gc:.2f}",
            },
            {
                "primer_name": "POOL_CONST_R",
                "sequence_5to3": reverse_primer,
                "orientation": "reverse",
                "adds_sequence": desired_3,
                "anneal_sequence": str(Seq(handle_3).reverse_complement()),
                "anneal_length_nt": len(handle_3),
                "anneal_tm_c": f"{right_tm:.2f}",
                "anneal_gc_percent": f"{right_gc:.2f}",
            },
        ]
        for primer in common_primers:
            ordering_rows.append(
                {
                    "name": primer["primer_name"],
                    "sequence": primer["sequence_5to3"],
                    "category": "common_primer",
                    "target_name": "ALL",
                    "order_scope": "order_once",
                }
            )

        pool_csv = output_dir / "synthetic_gene_pool.csv"
        primer_csv = output_dir / "synthetic_gene_pool_primers.csv"
        ordering_tsv = output_dir / "synthetic_gene_pool_ordering.tsv"
        instructions_txt = output_dir / "synthetic_gene_pool_instructions.txt"

        with pool_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "target_name",
                    "input_type",
                    "pool_oligo_name",
                    "sequence_5to3",
                    "pool_oligo_length",
                    "added_start_codon",
                    "added_stop_codon",
                    "warnings",
                ],
            )
            writer.writeheader()
            writer.writerows(pool_rows)

        with primer_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "primer_name",
                    "sequence_5to3",
                    "orientation",
                    "adds_sequence",
                    "anneal_sequence",
                    "anneal_length_nt",
                    "anneal_tm_c",
                    "anneal_gc_percent",
                ],
            )
            writer.writeheader()
            writer.writerows(common_primers)

        with ordering_tsv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["name", "sequence", "category", "target_name", "order_scope"])
            for row in ordering_rows:
                writer.writerow([row["name"], row["sequence"], row["category"], row["target_name"], row["order_scope"]])

        instructions_txt.write_text(
            "\n".join(
                [
                    "Synthetic Gene Pool Ordering Guidance",
                    "",
                    "This tool is specialized for E. coli T7 cell-free protein synthesis.",
                    "1. Order all pool oligos as the pooled synthesis set.",
                    "2. Order POOL_CONST_F and POOL_CONST_R once for the entire pool.",
                    "3. Amplify the pool with the common primer pair to lift out each gene while adding the final 5' and 3' sequences.",
                    "4. The forward primer adds the constant 5' region plus any N-terminal tag; the reverse primer adds any C-terminal tag plus stop codon and constant 3' region.",
                    f"5. Forward primer anneal handle length: {len(handle_5)} nt. Reverse primer anneal handle length: {len(handle_3)} nt.",
                    f"6. Common primer lengths: F={len(forward_primer)} nt, R={len(reverse_primer)} nt (must stay below 100 bp).",
                    "7. Ordered pool oligos are kept short by including only the shared anneal handles, not the full promoter/terminator payload.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        return {
            "pool_csv": pool_csv,
            "primer_csv": primer_csv,
            "ordering_tsv": ordering_tsv,
            "instructions_txt": instructions_txt,
            "run_log": run_log,
        }
    finally:
        if owns_logger:
            for handler in managed_logger.handlers[:]:
                handler.close()
                managed_logger.removeHandler(handler)

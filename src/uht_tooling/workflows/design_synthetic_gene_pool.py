"""Design pooled synthetic-gene oligos with reusable lift-out primers."""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from importlib import resources
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
    calc_tm,
)

CELL_FREE_TARGET_HOST = "e_coli"
MAX_COMMON_PRIMER_LENGTH = 99
PULLOUT_INDEX_LENGTH = 10


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


def _rewrite_pullout_three_prime_warnings(warnings: List[str]) -> List[str]:
    rewritten: List[str] = []
    for warning in warnings:
        if warning == "Added terminal stop codon on reusable 3' constant oligo.":
            rewritten.append("Added terminal stop codon on pooled ordering oligo upstream of the unique pullout index.")
        elif warning == "Moved terminal stop codon onto reusable 3' constant oligo.":
            rewritten.append("Placed terminal stop codon on pooled ordering oligo upstream of the unique pullout index.")
        else:
            rewritten.append(warning)
    return rewritten


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


def _build_pullout_three_prime_components(
    constant_3prime: str,
    stop_codon: str,
    tag_mode: str,
    tag_seq: str,
) -> Tuple[str, str, str, float, float]:
    three_prime_prefix = stop_codon
    if tag_mode.startswith("c_") and tag_seq:
        three_prime_prefix = tag_seq + three_prime_prefix
    handle_len, handle_tm, handle_gc = _pick_edge_overlap(constant_3prime, use_tail=False)
    terminator_handle = constant_3prime[:handle_len]
    return three_prime_prefix, constant_3prime, terminator_handle, handle_tm, handle_gc


@lru_cache(maxsize=1)
def _load_pullout_index_table() -> Dict[str, object]:
    table_path = resources.files("uht_tooling.data").joinpath("cfps_pullout_indexes.json")
    with table_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _assign_pullout_indexes(count: int) -> List[str]:
    table = _load_pullout_index_table()
    indexes = [str(entry["sequence"]) for entry in table["entries"]]
    if count > len(indexes):
        raise ValueError(
            f"Requested {count} pullout indexes, but only {len(indexes)} precomputed deterministic indexes are available in table version {table['table_version']}."
        )
    return indexes[:count]


def run_design_synthetic_gene_pool(
    sequence_fasta: Path,
    output_dir: Path,
    input_type: str = "auto",
    tag_mode: str = "none",
    include_pullout_primers: bool = False,
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
        pullout_index_table = _load_pullout_index_table() if include_pullout_primers else None
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
        pullout_three_prime_prefix = ""
        pullout_terminator = constant_3prime
        pullout_handle_3 = handle_3
        pullout_right_tm = right_tm
        pullout_right_gc = right_gc
        if include_pullout_primers:
            (
                pullout_three_prime_prefix,
                pullout_terminator,
                pullout_handle_3,
                pullout_right_tm,
                pullout_right_gc,
            ) = _build_pullout_three_prime_components(
                constant_3prime=constant_3prime,
                stop_codon=stop_codon,
                tag_mode=tag_mode,
                tag_seq=tag_seq,
            )

        pullout_indexes = _assign_pullout_indexes(len(records)) if include_pullout_primers else []

        pool_rows: List[Dict[str, object]] = []
        ordering_rows: List[Dict[str, object]] = []
        for idx, record in enumerate(records):
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
            if include_pullout_primers:
                warnings = _rewrite_pullout_three_prime_warnings(warnings)
            pool_seq = handle_5 + coding_seq + handle_3
            pool_name = f"{record.name}_POOL"
            warning_text = " | ".join(warnings)
            pullout_index = pullout_indexes[idx] if include_pullout_primers else ""
            if include_pullout_primers:
                pool_seq = handle_5 + coding_seq + pullout_three_prime_prefix + pullout_index + pullout_handle_3
            pool_rows.append(
                {
                    "target_name": record.name,
                    "input_type": resolved_input_type,
                    "pool_oligo_name": pool_name,
                    "sequence_5to3": pool_seq,
                    "pool_oligo_length": len(pool_seq),
                    "pullout_index": pullout_index,
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
        reverse_primer = str(Seq((pullout_terminator if include_pullout_primers else desired_3)).reverse_complement())
        if len(forward_primer) > MAX_COMMON_PRIMER_LENGTH or len(reverse_primer) > MAX_COMMON_PRIMER_LENGTH:
            raise ValueError(
                "The common primer pair exceeds the 99 bp ordering limit. Shorten the constant/tag payload for this cell-free design."
            )
        common_primers = [
            {
                "primer_name": "POOL_CONST_F",
                "sequence_5to3": forward_primer,
                "primer_role": "common_liftout",
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
                "primer_role": "common_liftout",
                "orientation": "reverse",
                "adds_sequence": pullout_terminator if include_pullout_primers else desired_3,
                "anneal_sequence": str(Seq((pullout_handle_3 if include_pullout_primers else handle_3)).reverse_complement()),
                "anneal_length_nt": len(pullout_handle_3 if include_pullout_primers else handle_3),
                "anneal_tm_c": f"{(pullout_right_tm if include_pullout_primers else right_tm):.2f}",
                "anneal_gc_percent": f"{(pullout_right_gc if include_pullout_primers else right_gc):.2f}",
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

        pullout_primers: List[Dict[str, object]] = []
        if include_pullout_primers:
            terminator_tail = pullout_terminator[len(pullout_handle_3) :]
            for row in pool_rows:
                pullout_index = str(row["pullout_index"])
                pullout_binding_region = pullout_index + pullout_handle_3
                pullout_primer = str(Seq(terminator_tail + pullout_binding_region).reverse_complement())
                if len(pullout_primer) > MAX_COMMON_PRIMER_LENGTH:
                    raise ValueError(
                        "The gene-specific pullout reverse primers exceed the 99 bp ordering limit. Shorten the constant/tag payload for this cell-free design."
                    )
                pullout_name = f"{row['target_name']}_PULLOUT_R"
                pullout_primers.append(
                    {
                        "primer_name": pullout_name,
                        "sequence_5to3": pullout_primer,
                        "primer_role": "gene_specific_pullout",
                        "orientation": "reverse",
                        "adds_sequence": terminator_tail,
                        "anneal_sequence": str(Seq(pullout_binding_region).reverse_complement()),
                        "anneal_length_nt": len(pullout_binding_region),
                        "anneal_tm_c": f"{calc_tm(pullout_binding_region):.2f}",
                        "anneal_gc_percent": f"{calc_gc_content(pullout_binding_region):.2f}",
                        "pullout_index": pullout_index,
                    }
                )
                ordering_rows.append(
                    {
                        "name": pullout_name,
                        "sequence": pullout_primer,
                        "category": "pullout_primer",
                        "target_name": row["target_name"],
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
                    "pullout_index",
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
                    "primer_role",
                    "orientation",
                    "adds_sequence",
                    "pullout_index",
                    "anneal_sequence",
                    "anneal_length_nt",
                    "anneal_tm_c",
                    "anneal_gc_percent",
                ],
            )
            writer.writeheader()
            writer.writerows(common_primers)
            writer.writerows(pullout_primers)

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
                    "3. Amplify the pool with the common primer pair to generate the full pooled CFPS library.",
                    "4. The forward primer adds the constant 5' region plus any N-terminal tag; the reverse primer completes the shared 3' payload for the pooled CFPS library.",
                    (
                        f"5. Forward primer anneal handle length: {len(handle_5)} nt. Reverse primer anneal handle length: "
                        f"{len(pullout_handle_3 if include_pullout_primers else handle_3)} nt."
                    ),
                    f"6. Common primer lengths: F={len(forward_primer)} nt, R={len(reverse_primer)} nt (must stay below 100 bp).",
                    "7. Ordered pool oligos are kept short by including only the shared anneal handles, not the full promoter/terminator payload.",
                ]
                + (
                    [
                        "8. Gene-specific pullout mode is enabled: each ordered pool oligo includes a deterministic unique index between the stop/c-tag segment and the terminator.",
                        "9. First amplify the full pooled library with POOL_CONST_F and POOL_CONST_R.",
                        "10. Then use POOL_CONST_F with the matching *_PULLOUT_R primer to selectively recover one variant from that full pooled CFPS library.",
                        (
                            "11. The gene-specific reverse primer binds the built-in "
                            f"{PULLOUT_INDEX_LENGTH} nt index plus the adjacent terminator handle; it does not introduce the unique index."
                        ),
                        f"12. Pullout index table version: {pullout_index_table['table_version']}.",
                    ]
                    if include_pullout_primers
                    else []
                )
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

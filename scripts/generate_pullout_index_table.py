"""Deterministically generate the packaged CFPS pullout-index table."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List

from Bio.Seq import Seq

DEFAULT_OUTPUT = Path("src/uht_tooling/data/cfps_pullout_indexes.json")
INDEX_LENGTH = 10
MIN_GC_COUNT = 4
MAX_GC_COUNT = 6
MAX_HOMOPOLYMER = 2
MIN_HAMMING_DISTANCE = 5
MAX_SUFFIX_PREFIX_OVERLAP = 2
MAX_SHARED_KMER = 5
MAX_SELF_SHARED_KMER = 3
PREFIX_BUCKET_LENGTH = 4


@dataclass(frozen=True)
class CandidateMetrics:
    sequence: str
    gc_percent: float
    reverse_complement: str
    shared_kmers_len6: frozenset[str]
    rc_kmers_len6: frozenset[str]


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def gc_percent(seq: str) -> float:
    return round(100.0 * sum(base in {"G", "C"} for base in seq) / len(seq), 2)


def hamming_distance(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def has_homopolymer(seq: str, max_run: int = MAX_HOMOPOLYMER) -> bool:
    run = 1
    for idx in range(1, len(seq)):
        run = run + 1 if seq[idx] == seq[idx - 1] else 1
        if run > max_run:
            return True
    return False


def kmers(seq: str, size: int) -> frozenset[str]:
    return frozenset(seq[idx : idx + size] for idx in range(len(seq) - size + 1))


def longest_suffix_prefix_overlap(left: str, right: str) -> int:
    max_len = min(len(left), len(right))
    for size in range(max_len, 0, -1):
        if left[-size:] == right[:size]:
            return size
    return 0


def longest_shared_kmer(left: str, right: str) -> int:
    max_len = min(len(left), len(right))
    for size in range(max_len, 0, -1):
        if kmers(left, size) & kmers(right, size):
            return size
    return 0


def enumerate_candidates() -> Iterable[CandidateMetrics]:
    for bases in product("ACGT", repeat=INDEX_LENGTH):
        sequence = "".join(bases)
        gc_count = sum(base in {"G", "C"} for base in sequence)
        if gc_count < MIN_GC_COUNT or gc_count > MAX_GC_COUNT:
            continue
        if has_homopolymer(sequence):
            continue
        rc_seq = reverse_complement(sequence)
        if longest_suffix_prefix_overlap(sequence, rc_seq) > MAX_SUFFIX_PREFIX_OVERLAP:
            continue
        if longest_suffix_prefix_overlap(rc_seq, sequence) > MAX_SUFFIX_PREFIX_OVERLAP:
            continue
        if longest_shared_kmer(sequence, rc_seq) > MAX_SELF_SHARED_KMER:
            continue
        yield CandidateMetrics(
            sequence=sequence,
            gc_percent=gc_percent(sequence),
            reverse_complement=rc_seq,
            shared_kmers_len6=kmers(sequence, MAX_SHARED_KMER + 1),
            rc_kmers_len6=kmers(rc_seq, MAX_SHARED_KMER + 1),
        )


def select_index_set(limit: int) -> List[CandidateMetrics]:
    buckets: Dict[str, List[CandidateMetrics]] = defaultdict(list)
    for candidate in enumerate_candidates():
        buckets[candidate.sequence[:PREFIX_BUCKET_LENGTH]].append(candidate)

    ordered_candidates: List[CandidateMetrics] = []
    bucket_keys = sorted(buckets)
    max_bucket_size = max((len(entries) for entries in buckets.values()), default=0)
    for round_idx in range(max_bucket_size):
        for bucket_key in bucket_keys:
            bucket = buckets[bucket_key]
            if round_idx < len(bucket):
                ordered_candidates.append(bucket[round_idx])

    selected: List[CandidateMetrics] = []
    for candidate in ordered_candidates:
        accepted = True
        for chosen in selected:
            if hamming_distance(candidate.sequence, chosen.sequence) < MIN_HAMMING_DISTANCE:
                accepted = False
                break
            if candidate.shared_kmers_len6 & chosen.shared_kmers_len6:
                accepted = False
                break
            if candidate.shared_kmers_len6 & chosen.rc_kmers_len6:
                accepted = False
                break
            if longest_suffix_prefix_overlap(candidate.sequence, chosen.sequence) > MAX_SUFFIX_PREFIX_OVERLAP:
                accepted = False
                break
            if longest_suffix_prefix_overlap(chosen.sequence, candidate.sequence) > MAX_SUFFIX_PREFIX_OVERLAP:
                accepted = False
                break
            if longest_suffix_prefix_overlap(candidate.sequence, chosen.reverse_complement) > MAX_SUFFIX_PREFIX_OVERLAP:
                accepted = False
                break
            if longest_suffix_prefix_overlap(chosen.reverse_complement, candidate.sequence) > MAX_SUFFIX_PREFIX_OVERLAP:
                accepted = False
                break
        if accepted:
            selected.append(candidate)
            if len(selected) >= limit:
                return selected
    return selected


def build_table(limit: int) -> Dict[str, object]:
    selected = select_index_set(limit)
    rows = []
    sequences = [entry.sequence for entry in selected]
    for entry in selected:
        other_sequences = [seq for seq in sequences if seq != entry.sequence]
        min_distance = min((hamming_distance(entry.sequence, seq) for seq in other_sequences), default=INDEX_LENGTH)
        rows.append(
            {
                "sequence": entry.sequence,
                "gc_percent": entry.gc_percent,
                "min_hamming_distance_to_table": min_distance,
                "screening_summary": {
                    "length_nt": INDEX_LENGTH,
                    "gc_count_range": [MIN_GC_COUNT, MAX_GC_COUNT],
                    "max_homopolymer_run": MAX_HOMOPOLYMER,
                    "max_self_shared_kmer_nt": MAX_SELF_SHARED_KMER,
                    "max_shared_kmer_with_table_nt": MAX_SHARED_KMER,
                    "max_suffix_prefix_overlap_nt": MAX_SUFFIX_PREFIX_OVERLAP,
                    "max_suffix_prefix_overlap_with_reverse_complement_nt": MAX_SUFFIX_PREFIX_OVERLAP,
                },
            }
        )
    return {
        "table_name": "cfps_pullout_indexes",
        "table_version": "2026-06-15-v2",
        "selection_method": "deterministic_greedy_prefix_bucket_round_robin",
        "selection_constraints": {
            "length_nt": INDEX_LENGTH,
            "gc_percent_range": [40, 60],
            "max_homopolymer_run": MAX_HOMOPOLYMER,
            "min_pairwise_hamming_distance": MIN_HAMMING_DISTANCE,
            "max_shared_kmer_with_table_nt": MAX_SHARED_KMER,
            "max_self_shared_kmer_nt": MAX_SELF_SHARED_KMER,
            "max_suffix_prefix_overlap_nt": MAX_SUFFIX_PREFIX_OVERLAP,
            "off_target_screening": "not_applied_context_independent_table",
            "candidate_order": f"prefix_round_robin_{PREFIX_BUCKET_LENGTH}nt_then_lexicographic_ACGT",
        },
        "entry_count": len(rows),
        "entries": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int, default=106)
    args = parser.parse_args()

    table = build_table(args.limit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(table, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

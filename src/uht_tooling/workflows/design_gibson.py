import argparse
import csv
import logging
import re
from collections import defaultdict
from math import floor
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

OVERHANG_LEN = 20
ANNEAL_LEN = 20


def codon_table() -> Dict[str, str]:
    return {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }


def translate_codon(cd: str) -> str:
    return codon_table().get(cd.upper(), "?")


def pick_mutant_codon(wt_codon: str, target_aa: str) -> str:
    best_list: List[tuple[str, int]] = []
    for codon, aa in codon_table().items():
        if aa == target_aa:
            diff = sum(a != b for a, b in zip(codon.upper(), wt_codon.upper()))
            best_list.append((codon.upper(), diff))
    if not best_list:
        return None
    best_list.sort(key=lambda x: x[1])
    return best_list[0][0]


def get_subseq_circ(seq: str, start: int, length: int) -> str:
    N = len(seq)
    if length >= N:
        raise ValueError(f"Requested length {length} ≥ sequence length {N}")
    s_mod = start % N
    end = s_mod + length
    if end <= N:
        return seq[s_mod:end]
    return seq[s_mod:] + seq[: (end % N)]


def design_gibson_primers(full_seq: str, region_start: int, old_len: int, new_seq: str):
    m = len(new_seq)
    if OVERHANG_LEN < m:
        raise ValueError(f"Mutation length {m} > OVERHANG_LEN ({OVERHANG_LEN})")

    N = len(full_seq)
    flank_left = floor((OVERHANG_LEN - m) / 2)
    flank_right = OVERHANG_LEN - m - flank_left

    if flank_left >= N or flank_right >= N or ANNEAL_LEN >= N:
        raise ValueError("Sequence too short for requested OVERHANG/ANNEAL lengths")

    oh_left = get_subseq_circ(full_seq, region_start - flank_left, flank_left)
    oh_right = get_subseq_circ(full_seq, region_start + old_len, flank_right)
    overhang = oh_left + new_seq + oh_right

    fwd_start = (region_start + old_len + flank_right) % N
    rev_start = (region_start - flank_left - ANNEAL_LEN) % N

    fwd_anneal = get_subseq_circ(full_seq, fwd_start, ANNEAL_LEN)
    rev_anneal = get_subseq_circ(full_seq, rev_start, ANNEAL_LEN)

    gibson_fwd = (overhang + fwd_anneal).upper()
    overhang_rc = str(Seq(overhang).reverse_complement()).upper()
    rev_anneal_rc = str(Seq(rev_anneal).reverse_complement()).upper()
    gibson_rev = (overhang_rc + rev_anneal_rc).upper()

    return gibson_fwd, gibson_rev, fwd_start, rev_start


def run_design_gibson(
    gene_fasta: Path,
    context_fasta: Path,
    mutations_csv: Path,
    output_dir: Path,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    gene_fasta = Path(gene_fasta)
    context_fasta = Path(context_fasta)
    mutations_csv = Path(mutations_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.design_gibson")
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
        gene_seq = str(next(SeqIO.parse(str(gene_fasta), "fasta")).seq).upper()
        context_seq = str(next(SeqIO.parse(str(context_fasta), "fasta")).seq).upper()
        logger.info("Loaded gene (%s nt) and context (%s nt).", len(gene_seq), len(context_seq))

        df = pd.read_csv(mutations_csv)
        if "mutations" not in df.columns:
            raise ValueError("Mutations CSV must contain a 'mutations' column.")
        entries = df["mutations"].dropna().tolist()
        logger.info("Loaded %s mutation entries.", len(entries))

        double_seq = context_seq + context_seq
        idx = double_seq.find(gene_seq)
        if idx == -1 or idx >= len(context_seq):
            raise ValueError("Could not align gene within circular context.")
        gene_offset = idx % len(context_seq)
        logger.info("Gene aligned at offset %s within context.", gene_offset)
        full_seq = context_seq

        primers_csv = output_dir / "Gibson_primers.csv"
        plan_csv = output_dir / "Gibson_assembly_plan.csv"
        group_entries: Dict[str, List[Dict[str, object]]] = defaultdict(list)

        with primers_csv.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Group", "Submutation", "Primer Name", "Sequence"])

            for entry in entries:
                submuts = entry.split("+")
                group_name = entry.replace("+", "_")
                logger.info("Processing group: %s with submutations: %s", group_name, submuts)

                for sub in submuts:
                    m_del = re.match(r"^([A-Z])(\d+)Del$", sub)
                    m_indel = re.match(r"^([A-Z])(\d+)InDel([A-Z])(\d+)([A-Z]+)$", sub)
                    m_sub = re.match(r"^([A-Z])(\d+)([A-Z])$", sub)
                    m_ins = re.match(r"^([A-Z])(\d+)([A-Z]{2,})$", sub)

                    if m_del:
                        _, pos1_s = m_del.groups()
                        region_start = gene_offset + (int(pos1_s) - 1) * 3
                        old_len = 3
                        new_seq = ""
                    elif m_indel:
                        wt1, pos1_s, _, pos2_s, ins_aa = m_indel.groups()
                        pos1, pos2 = int(pos1_s), int(pos2_s)
                        region_start = gene_offset + (pos1 - 1) * 3
                        old_len = (pos2 - pos1 + 1) * 3
                        wt_codon = get_subseq_circ(full_seq, region_start, 3)
                        new_seq = ""
                        for aa in ins_aa:
                            codon = pick_mutant_codon(wt_codon, aa)
                            if not codon:
                                raise ValueError(f"No codon found for {wt1}->{ins_aa}")
                            new_seq += codon
                    elif m_ins:
                        wt_aa, pos1_s, ins_str = m_ins.groups()
                        pos1 = int(pos1_s)
                        codon_start_old = gene_offset + (pos1 - 1) * 3
                        wt_codon = get_subseq_circ(full_seq, codon_start_old, 3)
                        if ins_str[0] == wt_aa:
                            inserted_aas = ins_str[1:]
                            region_start = codon_start_old + 3
                            old_len = 0
                        else:
                            inserted_aas = ins_str
                            region_start = codon_start_old
                            old_len = 3
                        new_seq = ""
                        for aa in inserted_aas:
                            codon = pick_mutant_codon(wt_codon, aa)
                            if not codon:
                                raise ValueError(f"No codon for insertion amino acid {aa}")
                            new_seq += codon
                    elif m_sub:
                        wt_aa, pos1_s, mut_aa = m_sub.groups()
                        pos1 = int(pos1_s)
                        region_start = gene_offset + (pos1 - 1) * 3
                        old_len = 3
                        wt_codon = get_subseq_circ(full_seq, region_start, 3)
                        translated = translate_codon(wt_codon)
                        if translated != wt_aa:
                            raise ValueError(
                                f"For {sub}: expected {wt_aa}, found {translated} at {wt_codon}"
                            )
                        new_seq = pick_mutant_codon(wt_codon, mut_aa)
                        if not new_seq:
                            raise ValueError(f"No minimal-change codon for {wt_aa}->{mut_aa}")
                    else:
                        raise ValueError(f"Unknown mutation format: {sub}")

                    gibson_fwd, gibson_rev, fwd_start, rev_start = design_gibson_primers(
                        full_seq, region_start, old_len, new_seq
                    )
                    primer_fwd_name = f"{group_name}__{sub}_Gibson_F"
                    primer_rev_name = f"{group_name}__{sub}_Gibson_R"

                    writer.writerow([group_name, sub, primer_fwd_name, gibson_fwd])
                    writer.writerow([group_name, sub, primer_rev_name, gibson_rev])

                    group_entries[group_name].append(
                        {
                            "sub": sub,
                            "fwd_name": primer_fwd_name,
                            "rev_name": primer_rev_name,
                            "fwd_pos": fwd_start % len(full_seq),
                            "rev_pos": rev_start % len(full_seq),
                            "fwd_seq": gibson_fwd,
                            "rev_seq": gibson_rev,
                        }
                    )

        with plan_csv.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(
                [
                    "Group",
                    "Submutation",
                    "PCR_Primer_Forward",
                    "PCR_Primer_Reverse",
                    "Tm (celsius)",
                    "Amplicon Size (bp)",
                ]
            )

            for group_name, entries in group_entries.items():
                sorted_forwards = sorted(entries, key=lambda e: e["fwd_pos"])
                sorted_reverses = sorted(entries, key=lambda e: e["rev_pos"])
                n = len(sorted_forwards)
                N = len(full_seq)
                for i in range(n):
                    f_entry = sorted_forwards[i]
                    r_entry = sorted_reverses[(i + 1) % n]

                    Tm_fwd = mt.Tm_NN(f_entry["fwd_seq"])
                    Tm_rev = mt.Tm_NN(r_entry["rev_seq"])
                    Tm_pair = min(Tm_fwd, Tm_rev)

                    fwd_start = f_entry["fwd_pos"]
                    rev_start = r_entry["rev_pos"]
                    rev_end = (rev_start + ANNEAL_LEN - 1) % N
                    if rev_end >= fwd_start:
                        amp_size = rev_end - fwd_start + 1
                    else:
                        amp_size = (N - fwd_start) + (rev_end + 1)

                    writer.writerow(
                        [
                            group_name,
                            f_entry["sub"],
                            f_entry["fwd_name"],
                            r_entry["rev_name"],
                            f"{Tm_pair:.1f}",
                            amp_size,
                        ]
                    )

        logger.info("Wrote Gibson outputs to %s and %s", primers_csv, plan_csv)
        return {"primers_csv": primers_csv, "plan_csv": plan_csv}
    finally:
        if managed_logger and logger:
            for handler in list(logger.handlers):
                handler.close()
                logger.removeHandler(handler)
            logger.propagate = True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Design Gibson assembly primers from user-provided inputs.")
    parser.add_argument("--gene-fasta", required=True, type=Path, help="Path to gene FASTA file.")
    parser.add_argument("--context-fasta", required=True, type=Path, help="Path to circular context FASTA file.")
    parser.add_argument(
        "--mutations-csv",
        required=True,
        type=Path,
        help="CSV with a 'mutations' column (use '+' to chain sub-mutations).",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to write result CSV files.")
    parser.add_argument("--log-path", default=None, type=Path, help="Optional log file path.")
    return parser


def main(argv: Optional[List[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_design_gibson(
        gene_fasta=args.gene_fasta,
        context_fasta=args.context_fasta,
        mutations_csv=args.mutations_csv,
        output_dir=args.output_dir,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

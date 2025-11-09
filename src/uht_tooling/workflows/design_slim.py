import argparse
import csv
import logging
import re
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

MIN_LEN = 20
MAX_LEN = 55
TARGET_TM = 60.0
MAX_TM = 70.0
UPSTREAM_15 = 12


def codon_table():
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


def translate_codon(cd):
    """Translate a 3-nt codon to its amino acid."""
    return codon_table().get(cd.upper(), "?")


def pick_mutant_codon(wt_codon, target_aa):
    best_list = []
    for codon, aa in codon_table().items():
        if aa == target_aa:
            diff = sum(a != b for a, b in zip(codon.upper(), wt_codon.upper()))
            best_list.append((codon.upper(), diff))
    if not best_list:
        return None
    best_list.sort(key=lambda x: x[1])
    return best_list[0][0]


def calc_tm(seq):
    return mt.Tm_NN(seq)


def pick_downstream_matched(dna_segment):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_segment) + 1):
        sub = dna_segment[:length]
        if len(sub) > MAX_LEN:
            break
        tm_value = calc_tm(sub)
        if tm_value <= MAX_TM:
            diff = abs(tm_value - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub


def pick_upstream_matched(dna_rev):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_rev) + 1):
        sub = dna_rev[:length]
        if len(sub) > MAX_LEN:
            break
        tm_value = calc_tm(sub)
        if tm_value <= MAX_TM:
            diff = abs(tm_value - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub


def design_long_forward(full_seq, region_start, old_len, new_codon):
    if region_start < UPSTREAM_15:
        raise ValueError("Not enough upstream for Lf (15 bp).")
    lf_up = full_seq[region_start - UPSTREAM_15 : region_start]
    down_seg = full_seq[region_start + old_len :]
    matched_down = pick_downstream_matched(down_seg)
    lf_seq = lf_up + new_codon + matched_down
    length = len(lf_seq)
    if length < MIN_LEN or length > MAX_LEN:
        raise ValueError(f"Lf length {length} not in [{MIN_LEN}..{MAX_LEN}].")
    lf_start = region_start - UPSTREAM_15
    lf_end = lf_start + length - 1
    return lf_seq.upper(), lf_start, lf_end


def design_long_reverse(full_seq, region_start, old_len, new_codon):
    up_rev = full_seq[:region_start][::-1]
    matched_up = pick_upstream_matched(up_rev)
    up_fwd = matched_up[::-1]
    lr_fwd = up_fwd + new_codon
    length = len(lr_fwd)
    if length < MIN_LEN or length > MAX_LEN:
        raise ValueError(f"Lr length {length} not in [{MIN_LEN}..{MAX_LEN}].")
    lr_start = region_start - len(up_fwd)
    lr_end = region_start + old_len - 1
    if lr_start < 0 or lr_end >= len(full_seq):
        raise ValueError("Lr coverage out of range.")
    lr_final = str(Seq(lr_fwd).reverse_complement()).upper()
    return lr_final, lr_start, lr_end


def design_short_reverse(full_seq, lf_start):
    sr_end = lf_start - 1
    if sr_end < 0:
        raise ValueError("No space for Sr (lf_start=0).")
    best_seg, best_diff = "", float("inf")
    for length in range(MIN_LEN, MAX_LEN + 1):
        start = sr_end - length + 1
        if start < 0:
            break
        region = full_seq[start : sr_end + 1]
        tm_value = calc_tm(region)
        if tm_value <= MAX_TM:
            diff = abs(tm_value - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, region
    if not best_seg:
        raise ValueError("Cannot find Sr [20..50] with Tm <=70 upstream.")
    return str(Seq(best_seg).reverse_complement()).upper()


def design_short_forward(full_seq, lr_end):
    sf_start = lr_end + 1
    if sf_start >= len(full_seq):
        raise ValueError("No space for Sf (lr_end near end).")
    best_seg, best_diff = "", float("inf")
    max_possible = len(full_seq) - sf_start
    for length in range(MIN_LEN, min(MAX_LEN, max_possible) + 1):
        sub = full_seq[sf_start : sf_start + length]
        tm_value = calc_tm(sub)
        if tm_value <= MAX_TM:
            diff = abs(tm_value - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, sub
    if not best_seg:
        raise ValueError("Cannot find Sf [20..50] with Tm <=70.")
    return best_seg.upper()


def run_design_slim(
    gene_fasta: Path,
    context_fasta: Path,
    mutations_csv: Path,
    output_dir: Path,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
):
    gene_fasta = Path(gene_fasta)
    context_fasta = Path(context_fasta)
    mutations_csv = Path(mutations_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.design_slim")
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
        gene_record = next(SeqIO.parse(str(gene_fasta), "fasta"))
        context_record = next(SeqIO.parse(str(context_fasta), "fasta"))
        gene = str(gene_record.seq).upper()
        context = str(context_record.seq).upper()
        logger.info("Loaded gene (%s nt) and context (%s nt).", len(gene), len(context))

        df = pd.read_csv(mutations_csv)
        if "mutations" not in df.columns:
            raise ValueError("Mutations CSV must contain a 'mutations' column.")
        mutations = df["mutations"].dropna().tolist()
        logger.info("Loaded %s mutation entries.", len(mutations))

        try:
            gene_offset = context.index(gene)
            logger.info("Gene aligned at offset %s within context.", gene_offset)
        except ValueError as exc:
            message = "Could not align gene within context. No perfect substring match found."
            logger.error(message)
            raise ValueError(message) from exc

        full_seq = context
        results_path = output_dir / "SLIM_primers.csv"
        with results_path.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Primer Name", "Sequence"])

            for mutation in mutations:
                try:
                    m = mutation
                    m_del = re.match(r"^([A-Z])(\d+)Del$", m)
                    m_indel = re.match(r"^([A-Z])(\d+)InDel([A-Z])(\d+)([A-Z]+)$", m)
                    m_sub = re.match(r"^([A-Z])(\d+)([A-Z])$", m)
                    m_ins = re.match(r"^([A-Z])(\d+)([A-Z]{2,})$", m)

                    if m_del:
                        wt_aa, pos1 = m_del.group(1), int(m_del.group(2))
                        region_start = gene_offset + (pos1 - 1) * 3
                        old_len = 3
                        new_seq = ""
                    elif m_indel:
                        wt1, pos1_s, wt2, pos2_s, ins_aa = m_indel.groups()
                        pos1, pos2 = int(pos1_s), int(pos2_s)
                        region_start = gene_offset + (pos1 - 1) * 3
                        old_len = (pos2 - pos1 + 1) * 3
                        wt_codon = full_seq[region_start : region_start + 3]
                        new_seq = ""
                        for aa in ins_aa:
                            codon = pick_mutant_codon(wt_codon, aa)
                            if not codon:
                                logger.error("No codon found for %s->%s", wt1, ins_aa)
                                raise ValueError(f"No codon found for {wt1}->{ins_aa}")
                            new_seq += codon
                    elif m_ins:
                        wt_aa, pos1_s, ins_str = m_ins.groups()
                        pos1 = int(pos1_s)
                        codon_start_old = gene_offset + (pos1 - 1) * 3
                        wt_codon = full_seq[codon_start_old : codon_start_old + 3]
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
                                logger.error("No codon for insertion amino acid %s", aa)
                                raise ValueError(f"No codon for insertion amino acid {aa}")
                            new_seq += codon
                    elif m_sub:
                        wt_aa, pos1, mut_aa = m_sub.group(1), int(m_sub.group(2)), m_sub.group(3)
                        region_start = gene_offset + (pos1 - 1) * 3
                        old_len = 3
                        wt_codon = full_seq[region_start : region_start + 3]
                        translated = translate_codon(wt_codon)
                        if translated != wt_aa:
                            logger.error(
                                "Expected %s but found %s at codon %s for mutation %s",
                                wt_aa,
                                translated,
                                wt_codon,
                                mutation,
                            )
                            raise ValueError(
                                "For {mutation}: expected {wt}, found {found} at {codon}".format(
                                    mutation=mutation, wt=wt_aa, found=translated, codon=wt_codon
                                )
                            )
                        new_seq = pick_mutant_codon(wt_codon, mut_aa)
                        if not new_seq:
                            logger.error("No minimal-change codon for %s->%s", wt_aa, mut_aa)
                            raise ValueError(f"No minimal-change codon for {wt_aa}->{mut_aa}")
                    else:
                        logger.error("Unknown mutation format: %s", mutation)
                        raise ValueError(f"Unknown mutation format: {mutation}")

                    lf_seq, lf_start, lf_end = design_long_forward(full_seq, region_start, old_len, new_seq)
                    lr_seq, lr_start, lr_end = design_long_reverse(full_seq, region_start, old_len, new_seq)
                    sr_seq = design_short_reverse(full_seq, lf_start)
                    sf_seq = design_short_forward(full_seq, lr_end)

                    writer.writerow([f"{mutation}_Lf", lf_seq])
                    writer.writerow([f"{mutation}_Sr", sr_seq])
                    writer.writerow([f"{mutation}_Lr", lr_seq])
                    writer.writerow([f"{mutation}_Sf", sf_seq])
                    logger.info("Designed primers for %s", mutation)
                except Exception as exc:
                    logger.error("Error processing mutation %s: %s", mutation, exc)
                    raise

        logger.info("SLIM primer design completed successfully. Output written to %s", results_path)
        return results_path
    finally:
        if managed_logger and logger:
            for handler in list(logger.handlers):
                handler.close()
                logger.removeHandler(handler)
            logger.propagate = True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Design SLIM primers from user-provided inputs.")
    parser.add_argument("--gene-fasta", required=True, type=Path, help="Path to FASTA file containing the gene sequence.")
    parser.add_argument(
        "--context-fasta",
        required=True,
        type=Path,
        help="Path to FASTA file containing the plasmid or genomic context.",
    )
    parser.add_argument(
        "--mutations-csv",
        required=True,
        type=Path,
        help="CSV file containing a 'mutations' column with each mutation specification.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where results and logs will be written.",
    )
    parser.add_argument(
        "--log-path",
        default=None,
        type=Path,
        help="Optional path for the run log (defaults to console logging).",
    )
    return parser


def main(argv: Optional[List[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_design_slim(
        gene_fasta=args.gene_fasta,
        context_fasta=args.context_fasta,
        mutations_csv=args.mutations_csv,
        output_dir=args.output_dir,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

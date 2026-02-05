"""KLD (Kinase-Ligation-DpnI) primer design for inverse PCR mutagenesis."""

import argparse
import csv
import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# KLD primer design parameters
MIN_TM = 50.0           # Minimum melting temperature
MAX_TM = 65.0           # Maximum melting temperature
TM_DIFF_THRESHOLD = 5.0  # Max Tm difference between F & R
MIN_LENGTH = 18         # Minimum primer length (template-binding region)
MAX_LENGTH = 24         # Maximum primer length (template-binding region)
MIN_GC = 40.0           # Minimum GC content %
MAX_GC = 60.0           # Maximum GC content %
MIN_BINDING = 10        # Minimum template-binding region

# IUPAC ambiguity codes mapping
IUPAC_AMBIGUITY = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'],      # puRine
    'Y': ['C', 'T'],      # pYrimidine
    'S': ['G', 'C'],      # Strong
    'W': ['A', 'T'],      # Weak
    'K': ['G', 'T'],      # Keto
    'M': ['A', 'C'],      # aMino
    'B': ['C', 'G', 'T'],  # not A
    'D': ['A', 'G', 'T'],  # not C
    'H': ['A', 'C', 'T'],  # not G
    'V': ['A', 'C', 'G'],  # not T
    'N': ['A', 'C', 'G', 'T'],
}

VALID_DEGENERATE_BASES = set(IUPAC_AMBIGUITY.keys())


def is_valid_degenerate_codon(codon: str) -> bool:
    """Check if a codon contains only valid IUPAC nucleotide codes."""
    return len(codon) == 3 and all(b.upper() in VALID_DEGENERATE_BASES for b in codon)


def contains_degenerate_bases(seq: str) -> bool:
    """Return True if sequence contains non-standard (degenerate) bases."""
    return any(b.upper() not in {'A', 'C', 'G', 'T'} for b in seq)


def expand_degenerate_sequence(seq: str) -> List[str]:
    """Expand a degenerate sequence to all possible standard sequences."""
    from itertools import product
    possibilities = [IUPAC_AMBIGUITY.get(b.upper(), [b]) for b in seq]
    return [''.join(combo) for combo in product(*possibilities)]


def codon_table():
    return {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }


def translate_codon(cd: str) -> str:
    """Translate a 3-nt codon to its amino acid."""
    return codon_table().get(cd.upper(), "?")


def pick_mutant_codon(wt_codon: str, target_aa: str) -> Optional[str]:
    """Pick the codon for target_aa that differs minimally from wt_codon."""
    best_list = []
    for codon, aa in codon_table().items():
        if aa == target_aa:
            diff = sum(a != b for a, b in zip(codon.upper(), wt_codon.upper()))
            best_list.append((codon.upper(), diff))
    if not best_list:
        return None
    best_list.sort(key=lambda x: x[1])
    return best_list[0][0]


def calc_tm_binding_region(seq: str) -> float:
    """Calculate Tm for template-binding region (handles degenerate bases)."""
    if contains_degenerate_bases(seq):
        expanded = expand_degenerate_sequence(seq)
        tms = [mt.Tm_NN(s) for s in expanded]
        return sum(tms) / len(tms)
    return mt.Tm_NN(seq)


def calc_gc_content(seq: str) -> float:
    """Calculate GC content percentage."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    return (gc / len(seq)) * 100 if seq else 0.0


def design_forward_primer(
    full_seq: str,
    mutation_nt_pos: int,
    new_codon: str,
    min_tm: float = MIN_TM,
    max_tm: float = MAX_TM,
    min_len: int = MIN_LENGTH,
    max_len: int = MAX_LENGTH,
) -> Tuple[str, str, float, float, int]:
    """
    Design forward primer with mutation at 5' end.

    Structure: [MUTATED_CODON] + [WT_DOWNSTREAM_SEQUENCE]

    The forward primer points downstream (→→→), with its 5' end at the mutation
    site. The template-binding region is the WT sequence immediately downstream
    of the mutation.

    Args:
        full_seq: Full plasmid/context sequence
        mutation_nt_pos: Nucleotide position of mutation start (0-indexed)
        new_codon: The mutated codon sequence (3 nt, may be degenerate)
        min_tm: Minimum melting temperature for binding region
        max_tm: Maximum melting temperature for binding region
        min_len: Minimum binding region length
        max_len: Maximum binding region length

    Returns:
        Tuple of (full_primer_seq, binding_region_seq, tm, gc, total_length)

    Raises:
        ValueError: If no valid binding region can be found
    """
    # Template-binding region starts immediately after mutation codon
    # (we use len(new_codon) to handle both regular codons and deletions)
    old_codon_len = 3  # Standard codon replacement
    binding_start = mutation_nt_pos + old_codon_len

    best_binding = None
    best_tm_diff = float('inf')
    target_tm = (min_tm + max_tm) / 2

    for length in range(min_len, max_len + 1):
        binding_end = binding_start + length
        if binding_end > len(full_seq):
            break

        binding_seq = full_seq[binding_start:binding_end]
        tm = calc_tm_binding_region(binding_seq)
        gc = calc_gc_content(binding_seq)

        if min_tm <= tm <= max_tm and MIN_GC <= gc <= MAX_GC:
            tm_diff = abs(tm - target_tm)
            if tm_diff < best_tm_diff:
                best_tm_diff = tm_diff
                best_binding = (binding_seq, tm, gc, length)

    if not best_binding:
        # Try to find any binding region that meets length constraints
        for length in range(min_len, max_len + 1):
            binding_end = binding_start + length
            if binding_end > len(full_seq):
                break
            binding_seq = full_seq[binding_start:binding_end]
            tm = calc_tm_binding_region(binding_seq)
            gc = calc_gc_content(binding_seq)
            # Relax constraints slightly
            if tm >= min_tm - 5 and tm <= max_tm + 5:
                best_binding = (binding_seq, tm, gc, length)
                break

    if not best_binding:
        raise ValueError(
            f"Cannot find forward primer binding region meeting constraints "
            f"(pos={mutation_nt_pos}, Tm={min_tm}-{max_tm}°C, len={min_len}-{max_len}bp)"
        )

    binding_seq, tm, gc, length = best_binding
    full_primer = new_codon + binding_seq
    return full_primer, binding_seq, tm, gc, len(full_primer)


def design_reverse_primer(
    full_seq: str,
    mutation_nt_pos: int,
    min_tm: float = MIN_TM,
    max_tm: float = MAX_TM,
    min_len: int = MIN_LENGTH,
    max_len: int = MAX_LENGTH,
) -> Tuple[str, float, float, int]:
    """
    Design reverse primer adjacent to forward primer's 5' end.

    Structure: reverse_complement([WT_UPSTREAM_SEQUENCE])

    The reverse primer's 5' end is at position mutation_nt_pos - 1, immediately
    adjacent to the forward primer's 5' end. The primer anneals to the top
    strand upstream of the mutation and points upstream (←←←).

    Args:
        full_seq: Full plasmid/context sequence
        mutation_nt_pos: Nucleotide position of mutation start (0-indexed)
        min_tm: Minimum melting temperature
        max_tm: Maximum melting temperature
        min_len: Minimum primer length
        max_len: Maximum primer length

    Returns:
        Tuple of (primer_seq, tm, gc, length)

    Raises:
        ValueError: If no valid primer can be found
    """
    # Upstream region ends at mutation position (exclusive)
    upstream_end = mutation_nt_pos

    best_primer = None
    best_tm_diff = float('inf')
    target_tm = (min_tm + max_tm) / 2

    for length in range(min_len, max_len + 1):
        upstream_start = upstream_end - length
        if upstream_start < 0:
            break

        upstream_seq = full_seq[upstream_start:upstream_end]
        tm = calc_tm_binding_region(upstream_seq)
        gc = calc_gc_content(upstream_seq)

        if min_tm <= tm <= max_tm and MIN_GC <= gc <= MAX_GC:
            tm_diff = abs(tm - target_tm)
            if tm_diff < best_tm_diff:
                best_tm_diff = tm_diff
                primer_seq = str(Seq(upstream_seq).reverse_complement())
                best_primer = (primer_seq, tm, gc, length)

    if not best_primer:
        # Try to find any primer that meets length constraints
        for length in range(min_len, max_len + 1):
            upstream_start = upstream_end - length
            if upstream_start < 0:
                break
            upstream_seq = full_seq[upstream_start:upstream_end]
            tm = calc_tm_binding_region(upstream_seq)
            gc = calc_gc_content(upstream_seq)
            # Relax constraints slightly
            if tm >= min_tm - 5 and tm <= max_tm + 5:
                primer_seq = str(Seq(upstream_seq).reverse_complement())
                best_primer = (primer_seq, tm, gc, length)
                break

    if not best_primer:
        raise ValueError(
            f"Cannot find reverse primer meeting constraints "
            f"(pos={mutation_nt_pos}, Tm={min_tm}-{max_tm}°C, len={min_len}-{max_len}bp)"
        )

    return best_primer


def balance_primer_tms(
    fwd_result: Tuple[str, str, float, float, int],
    rev_result: Tuple[str, float, float, int],
    full_seq: str,
    mutation_nt_pos: int,
    new_codon: str,
    tm_threshold: float = TM_DIFF_THRESHOLD,
) -> Tuple[Tuple[str, str, float, float, int], Tuple[str, float, float, int]]:
    """
    Balance Tm between forward and reverse primers by adjusting lengths.

    If the Tm difference exceeds the threshold, attempt to trim the hotter
    primer's binding region from its 3' end to reduce its Tm.

    Args:
        fwd_result: Forward primer tuple (full_seq, binding_seq, tm, gc, length)
        rev_result: Reverse primer tuple (primer_seq, tm, gc, length)
        full_seq: Full plasmid/context sequence
        mutation_nt_pos: Nucleotide position of mutation start
        new_codon: The mutated codon sequence
        tm_threshold: Maximum allowed Tm difference

    Returns:
        Adjusted (fwd_result, rev_result) tuples
    """
    fwd_seq, fwd_binding, fwd_tm, fwd_gc, fwd_len = fwd_result
    rev_seq, rev_tm, rev_gc, rev_len = rev_result

    tm_diff = abs(fwd_tm - rev_tm)

    if tm_diff <= tm_threshold:
        return fwd_result, rev_result

    # Try to trim the hotter primer
    if fwd_tm > rev_tm:
        # Trim forward binding region from 3' end
        while len(fwd_binding) > MIN_BINDING and fwd_tm - rev_tm > tm_threshold:
            fwd_binding = fwd_binding[:-1]
            fwd_tm = calc_tm_binding_region(fwd_binding)

        fwd_seq = new_codon + fwd_binding
        fwd_gc = calc_gc_content(fwd_seq)
        fwd_len = len(fwd_seq)
        fwd_result = (fwd_seq, fwd_binding, fwd_tm, fwd_gc, fwd_len)
    else:
        # Trim reverse primer from 3' end
        # The reverse primer is already reverse-complemented, so we need to
        # trim from the 3' end of the original upstream sequence
        upstream_end = mutation_nt_pos
        current_len = rev_len

        while current_len > MIN_BINDING and rev_tm - fwd_tm > tm_threshold:
            current_len -= 1
            upstream_start = upstream_end - current_len
            if upstream_start < 0:
                break
            upstream_seq = full_seq[upstream_start:upstream_end]
            rev_tm = calc_tm_binding_region(upstream_seq)

        if upstream_start >= 0:
            upstream_seq = full_seq[upstream_start:upstream_end]
            rev_seq = str(Seq(upstream_seq).reverse_complement())
            rev_gc = calc_gc_content(rev_seq)
            rev_len = len(rev_seq)
            rev_result = (rev_seq, rev_tm, rev_gc, rev_len)

    return fwd_result, rev_result


def run_design_kld(
    gene_fasta: Path,
    context_fasta: Path,
    mutations_csv: Path,
    output_dir: Path,
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Path:
    """
    Design KLD (inverse PCR) primers for mutations.

    KLD cloning uses two primers per mutation that point away from each other,
    amplifying the entire plasmid. The forward primer has the mutation at its
    5' end, and the reverse primer's 5' end is adjacent to the forward's.

    Args:
        gene_fasta: Path to FASTA file containing the gene sequence
        context_fasta: Path to FASTA file containing the plasmid/context sequence
        mutations_csv: CSV file with a 'mutations' column
        output_dir: Directory for output files
        log_path: Optional path for log file
        logger: Optional logger instance

    Returns:
        Path to the output CSV file

    Raises:
        ValueError: If inputs are invalid or mutations cannot be processed
    """
    gene_fasta = Path(gene_fasta)
    context_fasta = Path(context_fasta)
    mutations_csv = Path(mutations_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.design_kld")
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
        # Load sequences
        gene_record = next(SeqIO.parse(str(gene_fasta), "fasta"))
        context_record = next(SeqIO.parse(str(context_fasta), "fasta"))
        gene = str(gene_record.seq).upper()
        context = str(context_record.seq).upper()
        logger.info("Loaded gene (%s nt) and context (%s nt).", len(gene), len(context))

        # Load mutations
        df = pd.read_csv(mutations_csv)
        if "mutations" not in df.columns:
            raise ValueError("Mutations CSV must contain a 'mutations' column.")
        mutations = df["mutations"].dropna().tolist()
        logger.info("Loaded %s mutation entries.", len(mutations))

        # Align gene within context
        try:
            gene_offset = context.index(gene)
            logger.info("Gene aligned at offset %s within context.", gene_offset)
        except ValueError as exc:
            message = "Could not align gene within context. No perfect substring match found."
            logger.error(message)
            raise ValueError(message) from exc

        full_seq = context

        # Output file
        results_path = output_dir / "KLD_primers.csv"
        with results_path.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "Primer Name", "Sequence", "Tm (binding)", "GC%", "Length", "Notes"
            ])

            for mutation in mutations:
                try:
                    m = mutation
                    m_del = re.match(r"^([A-Z])(\d+)Del$", m)
                    m_indel = re.match(r"^([A-Z])(\d+)InDel([A-Z])(\d+)([A-Z]+)$", m)
                    m_sub = re.match(r"^([A-Z])(\d+)([A-Z])$", m)
                    m_ins = re.match(r"^([A-Z])(\d+)([A-Z]{2,})$", m)
                    m_lib = re.match(r"^([A-Z])(\d+):([A-Za-z]{3})$", m)

                    region_start: int
                    new_seq: str
                    notes: str = ""

                    if m_del:
                        # Deletion: remove the codon entirely
                        wt_aa, pos1 = m_del.group(1), int(m_del.group(2))
                        region_start = gene_offset + (pos1 - 1) * 3
                        # For deletion, new_seq is empty - forward primer starts at next codon
                        new_seq = ""
                        notes = "Deletion"

                    elif m_indel:
                        # InDel: replace range with new amino acids
                        wt1, pos1_s, wt2, pos2_s, ins_aa = m_indel.groups()
                        pos1, pos2 = int(pos1_s), int(pos2_s)
                        region_start = gene_offset + (pos1 - 1) * 3
                        wt_codon = full_seq[region_start : region_start + 3]
                        new_seq = ""
                        for aa in ins_aa:
                            codon = pick_mutant_codon(wt_codon, aa)
                            if not codon:
                                logger.error("No codon found for %s->%s", wt1, ins_aa)
                                raise ValueError(f"No codon found for {wt1}->{ins_aa}")
                            new_seq += codon
                        notes = f"InDel: {pos2 - pos1 + 1} AA -> {len(ins_aa)} AA"

                    elif m_ins:
                        # Insertion: add amino acids after position
                        wt_aa, pos1_s, ins_str = m_ins.groups()
                        pos1 = int(pos1_s)
                        codon_start_old = gene_offset + (pos1 - 1) * 3
                        wt_codon = full_seq[codon_start_old : codon_start_old + 3]
                        if ins_str[0] == wt_aa:
                            # First AA matches WT, so insert after
                            inserted_aas = ins_str[1:]
                            region_start = codon_start_old + 3
                        else:
                            # Replace WT with insertion
                            inserted_aas = ins_str
                            region_start = codon_start_old
                        new_seq = ""
                        for aa in inserted_aas:
                            codon = pick_mutant_codon(wt_codon, aa)
                            if not codon:
                                logger.error("No codon for insertion amino acid %s", aa)
                                raise ValueError(f"No codon for insertion amino acid {aa}")
                            new_seq += codon
                        notes = f"Insertion: +{len(inserted_aas)} AA"

                    elif m_sub:
                        # Substitution: single amino acid change
                        wt_aa, pos1, mut_aa = m_sub.group(1), int(m_sub.group(2)), m_sub.group(3)
                        region_start = gene_offset + (pos1 - 1) * 3
                        wt_codon = full_seq[region_start : region_start + 3]
                        translated = translate_codon(wt_codon)
                        if translated != wt_aa:
                            logger.error(
                                "Expected %s but found %s at codon %s for mutation %s",
                                wt_aa, translated, wt_codon, mutation,
                            )
                            raise ValueError(
                                f"For {mutation}: expected {wt_aa}, found {translated} at {wt_codon}"
                            )
                        new_seq = pick_mutant_codon(wt_codon, mut_aa)
                        if not new_seq:
                            logger.error("No minimal-change codon for %s->%s", wt_aa, mut_aa)
                            raise ValueError(f"No minimal-change codon for {wt_aa}->{mut_aa}")
                        notes = f"Substitution: {wt_aa}->{mut_aa}"

                    elif m_lib:
                        # Library mutation with degenerate codon
                        wt_aa, pos_str, degenerate_codon = m_lib.groups()
                        pos = int(pos_str)
                        degenerate_codon = degenerate_codon.upper()

                        if not is_valid_degenerate_codon(degenerate_codon):
                            raise ValueError(f"Invalid degenerate codon: {degenerate_codon}")

                        region_start = gene_offset + (pos - 1) * 3
                        wt_codon = full_seq[region_start : region_start + 3]
                        translated = translate_codon(wt_codon)
                        if translated != wt_aa:
                            logger.error(
                                "Expected %s but found %s at codon %s for mutation %s",
                                wt_aa, translated, wt_codon, mutation,
                            )
                            raise ValueError(
                                f"For {mutation}: expected {wt_aa}, found {translated} at {wt_codon}"
                            )

                        new_seq = degenerate_codon

                        # Log library coverage info
                        expanded_codons = expand_degenerate_sequence(degenerate_codon)
                        unique_aas = set(
                            translate_codon(c) for c in expanded_codons if translate_codon(c) != '?'
                        )
                        logger.info(
                            "Library mutation %s: %d possible codons, %d amino acids",
                            mutation, len(expanded_codons), len(unique_aas)
                        )
                        notes = f"Library: {len(expanded_codons)} codons, {len(unique_aas)} AAs"

                    else:
                        logger.error("Unknown mutation format: %s", mutation)
                        raise ValueError(f"Unknown mutation format: {mutation}")

                    # Handle deletion specially - forward primer starts at next position
                    if m_del:
                        # For deletion, the forward primer binding region starts
                        # immediately after the deleted codon
                        fwd_binding_start = region_start + 3

                        # Find optimal forward binding region
                        best_fwd = None
                        best_tm_diff = float('inf')
                        target_tm = (MIN_TM + MAX_TM) / 2

                        for length in range(MIN_LENGTH, MAX_LENGTH + 1):
                            binding_end = fwd_binding_start + length
                            if binding_end > len(full_seq):
                                break
                            binding_seq = full_seq[fwd_binding_start:binding_end]
                            tm = calc_tm_binding_region(binding_seq)
                            gc = calc_gc_content(binding_seq)
                            if MIN_TM <= tm <= MAX_TM and MIN_GC <= gc <= MAX_GC:
                                tm_diff = abs(tm - target_tm)
                                if tm_diff < best_tm_diff:
                                    best_tm_diff = tm_diff
                                    best_fwd = (binding_seq, binding_seq, tm, gc, length)

                        if not best_fwd:
                            raise ValueError(
                                f"Cannot find forward primer for deletion {mutation}"
                            )

                        fwd_seq, fwd_binding, fwd_tm, fwd_gc, fwd_len = best_fwd

                        # Reverse primer design is the same
                        rev_seq, rev_tm, rev_gc, rev_len = design_reverse_primer(
                            full_seq, region_start
                        )
                    else:
                        # Standard primer design
                        fwd_result = design_forward_primer(
                            full_seq, region_start, new_seq
                        )
                        rev_result = design_reverse_primer(full_seq, region_start)

                        # Balance Tms
                        fwd_result, rev_result = balance_primer_tms(
                            fwd_result, rev_result, full_seq, region_start, new_seq
                        )

                        fwd_seq, fwd_binding, fwd_tm, fwd_gc, fwd_len = fwd_result
                        rev_seq, rev_tm, rev_gc, rev_len = rev_result

                    # Write forward primer
                    writer.writerow([
                        f"{mutation}_F",
                        fwd_seq,
                        f"{fwd_tm:.1f}",
                        f"{fwd_gc:.1f}",
                        fwd_len,
                        notes,
                    ])

                    # Write reverse primer
                    writer.writerow([
                        f"{mutation}_R",
                        rev_seq,
                        f"{rev_tm:.1f}",
                        f"{rev_gc:.1f}",
                        rev_len,
                        "",
                    ])

                    logger.info(
                        "Designed KLD primers for %s: F_Tm=%.1f°C, R_Tm=%.1f°C, diff=%.1f°C",
                        mutation, fwd_tm, rev_tm, abs(fwd_tm - rev_tm)
                    )

                except Exception as exc:
                    logger.error("Error processing mutation %s: %s", mutation, exc)
                    raise

        logger.info("KLD primer design completed successfully. Output written to %s", results_path)
        return results_path

    finally:
        if managed_logger and logger:
            for handler in list(logger.handlers):
                handler.close()
                logger.removeHandler(handler)
            logger.propagate = True


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Design KLD (inverse PCR) primers from user-provided inputs."
    )
    parser.add_argument(
        "--gene-fasta",
        required=True,
        type=Path,
        help="Path to FASTA file containing the gene sequence.",
    )
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
    """Main entry point for command-line usage."""
    parser = build_parser()
    args = parser.parse_args(argv)
    run_design_kld(
        gene_fasta=args.gene_fasta,
        context_fasta=args.context_fasta,
        mutations_csv=args.mutations_csv,
        output_dir=args.output_dir,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()

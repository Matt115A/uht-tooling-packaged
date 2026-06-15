"""Design overlap-extension PCR oligos for IVTT-ready gene constructs."""

from __future__ import annotations

import csv
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from uht_tooling.config import load_config

DEFAULT_TARGET_OLIGO_LENGTH = 40
DEFAULT_MAX_ORDER_ONCE_LENGTH = 80
MIN_OLIGO_LENGTH = 20
MAX_OLIGO_LENGTH = 50
MIN_OVERLAP = 18
MAX_OVERLAP = 24
TARGET_OVERLAP_TM = 60.0
MIN_OVERLAP_TM = 58.0
MAX_OVERLAP_TM = 62.0
MIN_OVERLAP_GC = 35.0
MAX_OVERLAP_GC = 65.0
STANDARD_STOP_CODONS = {"TAA", "TAG", "TGA"}
DNA_BASES = {"A", "C", "G", "T"}
PROTEIN_BASES = set("ACDEFGHIKLMNPQRSTVWY*")
VALID_INPUT_TYPES = {"auto", "dna", "protein"}
VALID_TAG_MODES = {"none", "n_his", "c_his", "n_his_flag", "c_his_flag"}
CONSTANT_LABELS = {"constant_5prime", "constant_tag", "start", "stop", "constant_3prime"}
WELL_ROWS = "ABCDEFGH"
DEFAULT_TARGET_HOST = "e_coli"

# Built-in defaults for a generic E. coli T7 IVTT cassette.
DEFAULT_CONSTANT_5PRIME_DNA = "TAATACGACTCACTATAGGGAGAAGGAGATATACAT"
DEFAULT_CONSTANT_3PRIME_DNA = "GCTAGTTATTGCTCAGCGGT"
DEFAULT_STOP_CODON = "TAA"
DEFAULT_TAGS = {
    "n_his": "CATCATCATCATCATCAT",
    "c_his": "CATCATCATCATCATCAT",
    "n_his_flag": "CATCATCATCATCATCATGACTACAAAGACGATGACGACAAG",
    "c_his_flag": "CATCATCATCATCATCATGACTACAAAGACGATGACGACAAG",
}
DEFAULT_CODON_TABLE = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


@dataclass(frozen=True)
class SequencePart:
    name: str
    seq: str


@dataclass(frozen=True)
class Span:
    start: int
    end: int
    label: str


@dataclass(frozen=True)
class OligoTile:
    start: int
    end: int
    overlap_with_prev: int
    overlap_tm: Optional[float]
    overlap_gc: Optional[float]
    warning: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class TargetRecord:
    name: str
    raw_sequence: str


@dataclass(frozen=True)
class TargetDesign:
    target_name: str
    input_type: str
    construct_sequence: str
    tiles: Sequence[OligoTile]
    spans: Sequence[Span]
    warnings: Sequence[str]
    added_start: bool
    added_stop: bool
    tag_mode: str
    oligo_rows: Sequence[Dict[str, object]]
    assembly_strategy: str


def _setup_logger(log_path: Optional[Path], logger: Optional[logging.Logger]) -> Tuple[logging.Logger, bool]:
    if logger is not None:
        return logger, False
    managed_logger = logging.getLogger("uht_tooling.design_gene_oligos")
    managed_logger.setLevel(logging.INFO)
    handler: logging.Handler
    if log_path:
        handler = logging.FileHandler(log_path, mode="w")
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    managed_logger.handlers = []
    managed_logger.addHandler(handler)
    managed_logger.propagate = False
    return managed_logger, True


def _read_fasta_records(path: Path) -> List[TargetRecord]:
    records = [TargetRecord(record.id, str(record.seq)) for record in SeqIO.parse(str(path), "fasta")]
    if not records:
        raise ValueError("No FASTA records found in the input file.")
    return records


def _clean_sequence_text(seq: str) -> str:
    return "".join(seq.split()).upper()


def _detect_input_type(seq: str) -> str:
    letters = set(seq)
    if letters <= DNA_BASES:
        return "dna"
    if letters <= PROTEIN_BASES:
        return "protein"
    raise ValueError("Input sequence contains unsupported characters for DNA or protein input.")


def _normalize_dna(seq: str, logger: logging.Logger) -> Tuple[str, List[str], bool]:
    warnings: List[str] = []
    cleaned = _clean_sequence_text(seq)
    if not cleaned:
        raise ValueError("DNA sequence is empty.")
    invalid = set(cleaned) - DNA_BASES
    if invalid:
        raise ValueError(f"DNA input contains invalid bases: {''.join(sorted(invalid))}")
    if len(cleaned) % 3 != 0:
        raise ValueError("DNA input length must be divisible by 3.")
    had_start = cleaned.startswith("ATG")
    if cleaned[-3:] in STANDARD_STOP_CODONS:
        cleaned = cleaned[:-3]
        warnings.append("Stripped terminal stop codon from DNA input.")
        logger.warning(warnings[-1])
    if not cleaned:
        raise ValueError("DNA input only contained a stop codon.")
    protein = str(Seq(cleaned).translate())
    if "*" in protein:
        raise ValueError("DNA input contains an internal stop codon and does not look like a clean CDS.")
    return cleaned, warnings, had_start


def _normalize_protein(seq: str, logger: logging.Logger) -> Tuple[str, List[str], bool]:
    warnings: List[str] = []
    cleaned = _clean_sequence_text(seq)
    if not cleaned:
        raise ValueError("Protein sequence is empty.")
    invalid = set(cleaned) - PROTEIN_BASES
    if invalid:
        raise ValueError(f"Protein input contains invalid residues: {''.join(sorted(invalid))}")
    if cleaned.endswith("*"):
        cleaned = cleaned[:-1]
        warnings.append("Stripped trailing stop symbol from protein input.")
        logger.warning(warnings[-1])
    if "*" in cleaned:
        raise ValueError("Protein input contains an internal stop codon.")
    if not cleaned:
        raise ValueError("Protein input only contained a stop codon.")
    return cleaned, warnings, cleaned.startswith("M")


def _load_gene_oligo_config(config: Optional[Dict[str, object]]) -> Dict[str, object]:
    if config is None:
        config = load_config()
    gene_cfg = config.get("gene_oligos", {}) if isinstance(config, dict) else {}
    if not isinstance(gene_cfg, dict):
        raise ValueError("Config section 'gene_oligos' must be a mapping.")
    return gene_cfg


def _resolve_string(config: Dict[str, object], key: str, default: str) -> str:
    value = config.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Config value 'gene_oligos.{key}' must be a non-empty DNA string.")
    cleaned = _clean_sequence_text(value)
    invalid = set(cleaned) - DNA_BASES
    if invalid:
        raise ValueError(f"Config value 'gene_oligos.{key}' contains invalid DNA bases: {''.join(sorted(invalid))}")
    return cleaned


def _load_tag_sequence(config: Dict[str, object], tag_mode: str) -> str:
    if tag_mode == "none":
        return ""
    tags = config.get("tags", {})
    if tags and not isinstance(tags, dict):
        raise ValueError("Config section 'gene_oligos.tags' must be a mapping.")
    tag_entry = tags.get(tag_mode) if isinstance(tags, dict) else None
    if isinstance(tag_entry, dict):
        value = tag_entry.get("dna")
    else:
        value = DEFAULT_TAGS.get(tag_mode)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Tag sequence for '{tag_mode}' is not configured and has no built-in default.")
    cleaned = _clean_sequence_text(value)
    invalid = set(cleaned) - DNA_BASES
    if invalid:
        raise ValueError(
            f"Config value 'gene_oligos.tags.{tag_mode}.dna' contains invalid DNA bases: {''.join(sorted(invalid))}"
        )
    return cleaned


def _load_codon_table(config: Dict[str, object]) -> Dict[str, str]:
    codon_table = config.get("codon_table", DEFAULT_CODON_TABLE)
    if not isinstance(codon_table, dict) or not codon_table:
        raise ValueError("Config section 'gene_oligos.codon_table' must be a mapping.")
    normalized: Dict[str, str] = {}
    for aa, codon in codon_table.items():
        if not isinstance(aa, str) or not isinstance(codon, str):
            raise ValueError("Codon table entries must be string-to-string mappings.")
        aa_clean = aa.strip().upper()
        codon_clean = _clean_sequence_text(codon)
        if len(aa_clean) != 1 or aa_clean not in PROTEIN_BASES - {"*"}:
            raise ValueError(f"Invalid amino-acid key in codon table: {aa!r}")
        if len(codon_clean) != 3 or set(codon_clean) - DNA_BASES:
            raise ValueError(f"Invalid codon for amino acid {aa_clean}: {codon!r}")
        normalized[aa_clean] = codon_clean
    return normalized


def reverse_translate_protein(protein_seq: str, codon_table: Dict[str, str]) -> str:
    missing = sorted({aa for aa in protein_seq if aa not in codon_table})
    if missing:
        raise ValueError(f"Codon table is missing mappings for: {', '.join(missing)}")
    return "".join(codon_table[aa] for aa in protein_seq)


def codon_optimize_protein(protein_seq: str, target_host: str, codon_table: Dict[str, str]) -> str:
    try:
        from dnachisel import CodonOptimize, DnaOptimizationProblem, EnforceTranslation
    except ImportError as exc:  # pragma: no cover - exercised in runtime envs without dependency
        raise ValueError(
            "Protein input requires dnachisel for host-aware codon optimization. "
            "Install project dependencies to enable this workflow."
        ) from exc

    seed_sequence = reverse_translate_protein(protein_seq, codon_table)
    problem = DnaOptimizationProblem(
        sequence=seed_sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species=target_host)],
    )
    problem.resolve_constraints()
    problem.optimize()
    return str(problem.sequence)


def calc_gc_content(seq: str) -> float:
    gc = sum(1 for base in seq.upper() if base in {"G", "C"})
    return (gc / len(seq)) * 100 if seq else 0.0


def calc_tm(seq: str) -> float:
    return mt.Tm_NN(seq)


def _score_overlap(seq: str) -> Tuple[float, float, float, bool]:
    tm = calc_tm(seq)
    gc = calc_gc_content(seq)
    within = MIN_OVERLAP_TM <= tm <= MAX_OVERLAP_TM and MIN_OVERLAP_GC <= gc <= MAX_OVERLAP_GC
    penalty = abs(tm - TARGET_OVERLAP_TM)
    if gc < MIN_OVERLAP_GC:
        penalty += MIN_OVERLAP_GC - gc
    elif gc > MAX_OVERLAP_GC:
        penalty += gc - MAX_OVERLAP_GC
    return penalty, tm, gc, within


def _candidate_lengths(max_oligo_length: int) -> List[int]:
    return list(range(max_oligo_length, MIN_OLIGO_LENGTH - 1, -1))


def tile_sequence_with_overlaps(full_seq: str, oligo_length: int) -> List[OligoTile]:
    if not full_seq:
        raise ValueError("Final construct sequence is empty.")
    if len(full_seq) <= oligo_length:
        return [OligoTile(0, len(full_seq), 0, None, None, "")]

    memo: Dict[int, Optional[List[OligoTile]]] = {}

    def solve(start: int) -> Optional[List[OligoTile]]:
        if start in memo:
            return memo[start]
        remaining = len(full_seq) - start
        if MIN_OLIGO_LENGTH <= remaining <= oligo_length:
            memo[start] = [OligoTile(start, len(full_seq), 0, None, None, "")]
            return memo[start]

        end = start + oligo_length
        if end >= len(full_seq):
            memo[start] = None
            return None
        for overlap_len in range(MAX_OVERLAP, MIN_OVERLAP - 1, -1):
            overlap_start = end - overlap_len
            next_remaining = len(full_seq) - overlap_start
            if next_remaining < MIN_OLIGO_LENGTH:
                continue
            downstream = solve(overlap_start)
            if downstream is None:
                continue
            overlap_seq = full_seq[overlap_start:end]
            _, tm, gc, within = _score_overlap(overlap_seq)
            warning = ""
            if not within:
                warning = f"Overlap outside preferred window (Tm={tm:.1f}C, GC={gc:.1f}%)."
            tile = OligoTile(start, end, 0, None, None, warning)
            first = downstream[0]
            adjusted_first = OligoTile(
                first.start,
                first.end,
                overlap_len,
                tm,
                gc,
                first.warning,
            )
            memo[start] = [tile, adjusted_first, *downstream[1:]]
            return memo[start]

        memo[start] = None
        return None

    result = solve(0)
    if result is None:
        raise ValueError(
            f"Could not tile the construct into <= {oligo_length} nt oligos."
        )
    return result


def _build_sequence_parts(
    constant_5prime: str,
    tag_mode: str,
    tag_seq: str,
    goi_seq: str,
    stop_codon: str,
    constant_3prime: str,
) -> List[SequencePart]:
    parts: List[SequencePart] = [SequencePart("constant_5prime", constant_5prime)]
    parts.append(SequencePart("start", "ATG"))
    if tag_mode.startswith("n_") and tag_seq:
        parts.append(SequencePart("constant_tag", tag_seq))
    parts.append(SequencePart("goi", goi_seq))
    if tag_mode.startswith("c_") and tag_seq:
        parts.append(SequencePart("constant_tag", tag_seq))
    parts.append(SequencePart("stop", stop_codon))
    parts.append(SequencePart("constant_3prime", constant_3prime))
    return parts


def _annotate_spans(parts: Sequence[SequencePart]) -> List[Span]:
    spans: List[Span] = []
    cursor = 0
    for part in parts:
        if not part.seq:
            continue
        spans.append(Span(cursor, cursor + len(part.seq), part.name))
        cursor += len(part.seq)
    return spans


def _segment_type_for_interval(spans: Sequence[Span], start: int, end: int) -> str:
    labels = {span.label for span in spans if span.start < end and span.end > start}
    if len(labels) == 1:
        return next(iter(labels))
    return "mixed"


def _oligo_role(segment_type: str) -> Tuple[str, str]:
    if segment_type in CONSTANT_LABELS:
        return "constant", "order_once"
    return "gene_specific", "order_per_target"


def _oligo_sequence(full_seq: str, tile: OligoTile, index: int, start_with_forward: bool = True) -> Tuple[str, str]:
    seq = full_seq[tile.start:tile.end]
    forward = (index % 2 == 1) if start_with_forward else (index % 2 == 0)
    if forward:
        return seq, "forward"
    return str(Seq(seq).reverse_complement()), "reverse"


def _tile_middle_sequence(full_seq: str, middle_start: int, middle_end: int, max_oligo_length: int) -> List[OligoTile]:
    best_tiles: Optional[List[OligoTile]] = None
    for candidate in _candidate_lengths(max_oligo_length):
        tiles = tile_sequence_with_overlaps(full_seq[middle_start:middle_end], candidate)
        if len(tiles) % 2 == 0:
            return tiles
        if best_tiles is None:
            best_tiles = tiles
    if best_tiles is None:
        raise ValueError("Could not tile the gene-specific assembly region.")
    return best_tiles


def _design_external_primer(seq: str, direction: str) -> Tuple[str, float, float]:
    best: Optional[Tuple[float, str, float, float]] = None
    for length in range(min(MAX_OVERLAP, len(seq)), MIN_OLIGO_LENGTH - 1, -1):
        if direction == "forward":
            candidate = seq[:length]
        else:
            candidate = str(Seq(seq[-length:]).reverse_complement())
        tm = calc_tm(candidate)
        gc = calc_gc_content(candidate)
        penalty = abs(tm - TARGET_OVERLAP_TM)
        if best is None or penalty < best[0]:
            best = (penalty, candidate, tm, gc)
    if best is None:
        raise ValueError("Constant cassette is too short to design reusable external primers.")
    return best[1], best[2], best[3]


def _pick_edge_overlap(seq: str, use_tail: bool) -> Tuple[int, float, float]:
    best: Optional[Tuple[float, int, float, float]] = None
    for length in range(MAX_OVERLAP, MIN_OVERLAP - 1, -1):
        candidate = seq[-length:] if use_tail else seq[:length]
        penalty, tm, gc, within = _score_overlap(candidate)
        if best is None or (within and not (MIN_OVERLAP_TM <= best[2] <= MAX_OVERLAP_TM and MIN_OVERLAP_GC <= best[3] <= MAX_OVERLAP_GC)) or penalty < best[0]:
            best = (penalty, length, tm, gc)
    if best is None:
        raise ValueError("Could not determine edge overlap.")
    return best[1], best[2], best[3]


def _constant_edge_sequences(tag_mode: str, constant_5prime: str, tag_seq: str, stop_codon: str, constant_3prime: str) -> Tuple[str, str]:
    left = constant_5prime + "ATG"
    right = stop_codon + constant_3prime
    if tag_mode.startswith("n_") and tag_seq:
        left += tag_seq
    if tag_mode.startswith("c_") and tag_seq:
        right = tag_seq + right
    return left, right


def _write_fasta(path: Path, entries: Iterable[Tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in entries:
            handle.write(f">{name}\n")
            for idx in range(0, len(seq), 80):
                handle.write(seq[idx: idx + 80] + "\n")


def _well_name(index: int) -> Tuple[str, str]:
    plate = index // 96 + 1
    pos = index % 96
    row = WELL_ROWS[pos % 8]
    col = pos // 8 + 1
    return f"Plate{plate}", f"{row}{col}"


def _assembly_strategy(oligo_count: int) -> str:
    if oligo_count <= 8:
        return "single_pool_oe_pcr"
    if oligo_count <= 16:
        return "two_block_preassembly_then_full_length_pcr"
    return "three_block_preassembly_then_full_length_pcr"


def _normalize_target(
    record: TargetRecord,
    resolved_input_type: str,
    tag_mode: str,
    tag_seq: str,
    stop_codon: str,
    codon_table: Dict[str, str],
    target_host: str,
    logger: logging.Logger,
) -> Tuple[str, List[str], bool, bool]:
    warnings: List[str] = []
    added_start = False
    added_stop = False

    if resolved_input_type == "dna":
        coding_seq, local_warnings, has_start = _normalize_dna(record.raw_sequence, logger)
        warnings.extend(local_warnings)
        had_stop = any("Stripped terminal stop codon" in warning for warning in warnings)
        if has_start:
            coding_seq = coding_seq[3:]
            warnings.append("Moved 5' start codon onto reusable 5' constant oligo.")
            logger.warning("%s: %s", record.name, warnings[-1])
        else:
            warnings.append("Added 5' start codon on reusable 5' constant oligo.")
            logger.warning("%s: %s", record.name, warnings[-1])
            added_start = True
        added_stop = not had_stop
        if added_stop:
            warnings.append("Added terminal stop codon on reusable 3' constant oligo.")
            logger.warning("%s: %s", record.name, warnings[-1])
        else:
            warnings.append("Moved terminal stop codon onto reusable 3' constant oligo.")
            logger.warning("%s: %s", record.name, warnings[-1])
        return coding_seq, warnings, added_start, added_stop

    protein_seq, local_warnings, has_start = _normalize_protein(record.raw_sequence, logger)
    warnings.extend(local_warnings)
    had_stop = any("Stripped trailing stop symbol" in warning for warning in warnings)
    if has_start:
        protein_seq = protein_seq[1:]
        warnings.append("Moved initiating methionine onto reusable 5' constant oligo.")
        logger.warning("%s: %s", record.name, warnings[-1])
    else:
        warnings.append("Added initiating methionine on reusable 5' constant oligo.")
        logger.warning("%s: %s", record.name, warnings[-1])
        added_start = True
    added_stop = not had_stop
    if added_stop:
        warnings.append("Added terminal stop codon on reusable 3' constant oligo.")
        logger.warning("%s: %s", record.name, warnings[-1])
    else:
        warnings.append("Moved terminal stop codon onto reusable 3' constant oligo.")
        logger.warning("%s: %s", record.name, warnings[-1])
    optimized = codon_optimize_protein(protein_seq, target_host=target_host, codon_table=codon_table)
    warnings.append(f"Applied host-aware codon optimization for protein input using target host '{target_host}'.")
    logger.info("%s: %s", record.name, warnings[-1])
    return optimized, warnings, added_start, added_stop


def _design_target(
    record: TargetRecord,
    resolved_input_type: str,
    constant_5prime: str,
    constant_3prime: str,
    left_constant_seq: str,
    right_constant_seq: str,
    stop_codon: str,
    tag_mode: str,
    tag_seq: str,
    codon_table: Dict[str, str],
    target_host: str,
    target_oligo_length: int,
    logger: logging.Logger,
) -> TargetDesign:
    goi_seq, warnings, added_start, added_stop = _normalize_target(
        record=record,
        resolved_input_type=resolved_input_type,
        tag_mode=tag_mode,
        tag_seq=tag_seq,
        stop_codon=stop_codon,
        codon_table=codon_table,
        target_host=target_host,
        logger=logger,
    )
    parts = _build_sequence_parts(
        constant_5prime=constant_5prime,
        tag_mode=tag_mode,
        tag_seq=tag_seq,
        goi_seq=goi_seq,
        stop_codon=stop_codon,
        constant_3prime=constant_3prime,
    )
    spans = _annotate_spans(parts)
    construct = "".join(part.seq for part in parts)
    left_overlap_len, left_overlap_tm, left_overlap_gc = _pick_edge_overlap(left_constant_seq, use_tail=True)
    right_overlap_len, right_overlap_tm, right_overlap_gc = _pick_edge_overlap(right_constant_seq, use_tail=False)
    middle_start = len(left_constant_seq) - left_overlap_len
    middle_end = len(construct) - len(right_constant_seq) + right_overlap_len
    if middle_end <= middle_start:
        raise ValueError("Constant edge oligos consume the entire construct; no gene-specific assembly region remains.")
    middle_tiles = _tile_middle_sequence(construct, middle_start, middle_end, target_oligo_length)
    tiles: List[OligoTile] = []
    for idx, tile in enumerate(middle_tiles):
        warning = tile.warning
        overlap_with_prev = tile.overlap_with_prev
        overlap_tm = tile.overlap_tm
        overlap_gc = tile.overlap_gc
        if idx == 0:
            overlap_with_prev = left_overlap_len
            overlap_tm = left_overlap_tm
            overlap_gc = left_overlap_gc
        absolute_tile = OligoTile(
            middle_start + tile.start,
            middle_start + tile.end,
            overlap_with_prev,
            overlap_tm,
            overlap_gc,
            warning,
        )
        tiles.append(absolute_tile)
    if tiles:
        last = tiles[-1]
        last_warning = last.warning
        tiles[-1] = OligoTile(
            last.start,
            last.end,
            last.overlap_with_prev,
            right_overlap_tm,
            right_overlap_gc,
            last_warning,
        )
    rows: List[Dict[str, object]] = []
    start_with_forward = len(tiles) % 2 == 1
    for idx, tile in enumerate(tiles, start=1):
        oligo_seq, orientation = _oligo_sequence(construct, tile, idx, start_with_forward=start_with_forward)
        segment_type = _segment_type_for_interval(spans, tile.start, tile.end)
        role, order_scope = ("gene_specific", "order_per_target")
        row_warnings = []
        if idx == 1:
            row_warnings.extend(warnings)
            row_warnings.append("Overlaps reusable 5' constant edge oligo.")
            if not start_with_forward:
                row_warnings.append("Starts in reverse orientation to alternate with the reusable 5' edge oligo.")
        if tile.warning:
            row_warnings.append(tile.warning)
        if idx == len(tiles):
            row_warnings.append("Overlaps reusable 3' constant edge oligo.")
        rows.append(
            {
                "target_name": record.name,
                "oligo_name": f"{record.name}_OE_{idx:02d}",
                "sequence_5to3": oligo_seq,
                "orientation": orientation,
                "start_1based": tile.start + 1,
                "end_1based": tile.end,
                "oligo_length": tile.length,
                "overlap_with_prev_nt": tile.overlap_with_prev,
                "overlap_tm_c": "" if tile.overlap_tm is None else f"{tile.overlap_tm:.2f}",
                "overlap_gc_percent": "" if tile.overlap_gc is None else f"{tile.overlap_gc:.2f}",
                "segment_type": segment_type,
                "oligo_role": role,
                "order_scope": order_scope,
                "warnings": " | ".join(row_warnings),
            }
        )
    return TargetDesign(
        target_name=record.name,
        input_type=resolved_input_type,
        construct_sequence=construct,
        tiles=tiles,
        spans=spans,
        warnings=warnings,
        added_start=added_start,
        added_stop=added_stop,
        tag_mode=tag_mode,
        oligo_rows=rows,
        assembly_strategy=_assembly_strategy(len(tiles)),
    )


def _render_assembly_instructions(
    target_designs: Sequence[TargetDesign],
    constant_oligo_count: int,
    external_forward_name: str,
    external_reverse_name: str,
) -> str:
    lines = [
        "Gene Oligo Assembly Guidance",
        "",
        "Ordering and stock setup",
        "- Order oligos in 96-well format at 20 uM.",
        "- Constant oligos and external primers are marked as order_once and only need to be ordered once.",
        "- Gene-specific oligos are marked as order_per_target and should be ordered for each target construct.",
        "- Prepare 1:1 dilutions to 10 uM working stocks before assembly.",
        "",
        "Recommended assembly workflow",
        f"- Reusable constant oligos in this campaign: {constant_oligo_count}.",
        f"- Use external primers {external_forward_name} and {external_reverse_name} for production-scale amplification of the final construct.",
        "- First-stage overlap-extension PCR should use equimolar oligos at roughly 10-25 nM final concentration per oligo.",
        "- For single-pool assemblies, combine all OE oligos in one reaction without external primers for the first cycles, then spike in the external primer pair for full-length amplification or re-amplify in a second PCR.",
        "- For block assemblies, preassemble each block separately, gel-check or size-check the blocks, then combine equimolar purified blocks with the external primer pair for final amplification.",
        "",
        "Per-target strategy",
    ]
    for design in target_designs:
        lines.append(
            f"- {design.target_name}: {len(design.tiles)} OE oligos, strategy={design.assembly_strategy}."
        )
    return "\n".join(lines) + "\n"


def run_design_gene_oligos(
    sequence_fasta: Path,
    output_dir: Path,
    input_type: str = "auto",
    target_oligo_length: int = DEFAULT_TARGET_OLIGO_LENGTH,
    max_order_once_length: int = DEFAULT_MAX_ORDER_ONCE_LENGTH,
    target_host: str = DEFAULT_TARGET_HOST,
    tag_mode: str = "none",
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    config: Optional[Dict[str, object]] = None,
) -> Dict[str, Path]:
    if input_type not in VALID_INPUT_TYPES:
        raise ValueError(f"input_type must be one of {sorted(VALID_INPUT_TYPES)}")
    if tag_mode not in VALID_TAG_MODES:
        raise ValueError(f"tag_mode must be one of {sorted(VALID_TAG_MODES)}")
    if target_oligo_length < MIN_OLIGO_LENGTH:
        raise ValueError(f"target_oligo_length must be at least {MIN_OLIGO_LENGTH}")
    if max_order_once_length < MIN_OLIGO_LENGTH:
        raise ValueError(f"max_order_once_length must be at least {MIN_OLIGO_LENGTH}")

    sequence_fasta = Path(sequence_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_log = Path(log_path) if log_path else output_dir / "run.log"
    managed_logger, owns_logger = _setup_logger(run_log, logger)

    try:
        records = _read_fasta_records(sequence_fasta)
        gene_cfg = _load_gene_oligo_config(config)
        constant_5prime = _resolve_string(gene_cfg, "constant_5prime_dna", DEFAULT_CONSTANT_5PRIME_DNA)
        constant_3prime = _resolve_string(gene_cfg, "constant_3prime_dna", DEFAULT_CONSTANT_3PRIME_DNA)
        stop_codon = _resolve_string(gene_cfg, "stop_codon", DEFAULT_STOP_CODON)
        if len(stop_codon) != 3 or stop_codon not in STANDARD_STOP_CODONS:
            raise ValueError("Config value 'gene_oligos.stop_codon' must be one of TAA, TAG, or TGA.")
        tag_seq = _load_tag_sequence(gene_cfg, tag_mode)
        codon_table = _load_codon_table(gene_cfg)
        left_constant_seq, right_constant_seq = _constant_edge_sequences(
            tag_mode=tag_mode,
            constant_5prime=constant_5prime,
            tag_seq=tag_seq,
            stop_codon=stop_codon,
            constant_3prime=constant_3prime,
        )
        if len(left_constant_seq) > max_order_once_length or len(right_constant_seq) > max_order_once_length:
            raise ValueError(
                "Reusable constant oligos exceed max_order_once_length. Increase the hidden extra setting or shorten the constant cassette."
            )

        target_designs: List[TargetDesign] = []
        for record in records:
            cleaned = _clean_sequence_text(record.raw_sequence)
            resolved_input_type = _detect_input_type(cleaned) if input_type == "auto" else input_type
            managed_logger.info("Designing target %s as %s input", record.name, resolved_input_type)
            target_designs.append(
                _design_target(
                    record=record,
                    resolved_input_type=resolved_input_type,
                    constant_5prime=constant_5prime,
                    constant_3prime=constant_3prime,
                    left_constant_seq=left_constant_seq,
                    right_constant_seq=right_constant_seq,
                    stop_codon=stop_codon,
                    tag_mode=tag_mode,
                    tag_seq=tag_seq,
                    codon_table=codon_table,
                    target_host=target_host,
                    target_oligo_length=target_oligo_length,
                    logger=managed_logger,
                )
            )

        external_forward_seq, external_forward_tm, external_forward_gc = _design_external_primer(
            constant_5prime, "forward"
        )
        external_reverse_seq, external_reverse_tm, external_reverse_gc = _design_external_primer(
            constant_3prime, "reverse"
        )

        oligo_csv = output_dir / "gene_oligos.csv"
        assembly_csv = output_dir / "assembly_report.csv"
        constructs_fasta = output_dir / "final_constructs.fasta"
        order_fasta = output_dir / "ordered_oligos.fasta"
        plate_csv = output_dir / "ordering_plate.csv"
        external_csv = output_dir / "external_primers.csv"
        instructions_txt = output_dir / "assembly_instructions.txt"

        all_rows: List[Dict[str, object]] = []
        constant_order_items: Dict[Tuple[str, str], Dict[str, object]] = {}
        gene_specific_order_items: List[Dict[str, object]] = []

        constant_edge_rows = [
            {
                "target_name": "ALL",
                "oligo_name": "CONST_EDGE_5P",
                "sequence_5to3": left_constant_seq,
                "orientation": "forward",
                "start_1based": "",
                "end_1based": "",
                "oligo_length": len(left_constant_seq),
                "overlap_with_prev_nt": "",
                "overlap_tm_c": "",
                "overlap_gc_percent": "",
                "segment_type": "constant_edge_5prime",
                "oligo_role": "constant",
                "order_scope": "order_once",
                "warnings": "Reusable 5' constant edge oligo carrying promoter/start and any N-terminal tag.",
            },
            {
                "target_name": "ALL",
                "oligo_name": "CONST_EDGE_3P",
                "sequence_5to3": str(Seq(right_constant_seq).reverse_complement()),
                "orientation": "reverse",
                "start_1based": "",
                "end_1based": "",
                "oligo_length": len(right_constant_seq),
                "overlap_with_prev_nt": "",
                "overlap_tm_c": "",
                "overlap_gc_percent": "",
                "segment_type": "constant_edge_3prime",
                "oligo_role": "constant",
                "order_scope": "order_once",
                "warnings": "Reusable 3' constant edge oligo carrying any C-terminal tag, stop codon, and terminator.",
            },
        ]
        all_rows.extend(constant_edge_rows)
        constant_order_items[(left_constant_seq, "forward")] = {
            "oligo_name": "CONST_EDGE_5P",
            "sequence_5to3": left_constant_seq,
            "orientation": "forward",
            "category": "constant",
            "target_name": "ALL",
        }
        constant_order_items[(str(Seq(right_constant_seq).reverse_complement()), "reverse")] = {
            "oligo_name": "CONST_EDGE_3P",
            "sequence_5to3": str(Seq(right_constant_seq).reverse_complement()),
            "orientation": "reverse",
            "category": "constant",
            "target_name": "ALL",
        }

        for design in target_designs:
            for row in design.oligo_rows:
                all_rows.append(row)
                key = (str(row["sequence_5to3"]), str(row["orientation"]))
                if row["oligo_role"] == "constant":
                    constant_order_items.setdefault(
                        key,
                        {
                            "oligo_name": f"CONST_OE_{len(constant_order_items) + 1:02d}",
                            "sequence_5to3": row["sequence_5to3"],
                            "orientation": row["orientation"],
                            "category": "constant",
                            "target_name": "ALL",
                        },
                    )
                else:
                    gene_specific_order_items.append(
                        {
                            "oligo_name": row["oligo_name"],
                            "sequence_5to3": row["sequence_5to3"],
                            "orientation": row["orientation"],
                            "category": "gene_specific",
                            "target_name": row["target_name"],
                        }
                    )

        constant_name_map = {
            (item["sequence_5to3"], item["orientation"]): item["oligo_name"]
            for item in constant_order_items.values()
        }
        for row in all_rows:
            if row["oligo_role"] == "constant":
                row["oligo_name"] = constant_name_map[(row["sequence_5to3"], row["orientation"])]

        with assembly_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "target_name",
                    "input_type",
                    "tag_mode",
                    "construct_length_nt",
                    "oligo_count",
                    "target_oligo_length_nt",
                    "added_start_codon",
                    "added_stop_codon",
                    "assembly_strategy",
                ]
            )
            for design in target_designs:
                writer.writerow(
                    [
                        design.target_name,
                        design.input_type,
                        design.tag_mode,
                        len(design.construct_sequence),
                        len(design.tiles),
                        target_oligo_length,
                        "yes" if design.added_start else "no",
                        "yes" if design.added_stop else "no",
                        design.assembly_strategy,
                    ]
                )

        _write_fasta(
            constructs_fasta,
            ((design.target_name, design.construct_sequence) for design in target_designs),
        )

        external_rows = [
            {
                "oligo_name": "CONST_EXT_F",
                "sequence_5to3": external_forward_seq,
                "orientation": "forward",
                "category": "external",
                "target_name": "ALL",
                "tm_c": f"{external_forward_tm:.2f}",
                "gc_percent": f"{external_forward_gc:.2f}",
            },
            {
                "oligo_name": "CONST_EXT_R",
                "sequence_5to3": external_reverse_seq,
                "orientation": "reverse",
                "category": "external",
                "target_name": "ALL",
                "tm_c": f"{external_reverse_tm:.2f}",
                "gc_percent": f"{external_reverse_gc:.2f}",
            },
        ]

        for external_row in external_rows:
            all_rows.append(
                {
                    "target_name": external_row["target_name"],
                    "oligo_name": external_row["oligo_name"],
                    "sequence_5to3": external_row["sequence_5to3"],
                    "orientation": external_row["orientation"],
                    "start_1based": "",
                    "end_1based": "",
                    "oligo_length": len(str(external_row["sequence_5to3"])),
                    "overlap_with_prev_nt": "",
                    "overlap_tm_c": external_row["tm_c"],
                    "overlap_gc_percent": external_row["gc_percent"],
                    "segment_type": "external_primer",
                    "oligo_role": "constant",
                    "order_scope": "order_once",
                    "warnings": "Reusable external amplification primer.",
                }
            )

        with oligo_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "target_name",
                    "oligo_name",
                    "sequence_5to3",
                    "orientation",
                    "start_1based",
                    "end_1based",
                    "oligo_length",
                    "overlap_with_prev_nt",
                    "overlap_tm_c",
                    "overlap_gc_percent",
                    "segment_type",
                    "oligo_role",
                    "order_scope",
                    "warnings",
                ],
            )
            writer.writeheader()
            writer.writerows(all_rows)

        with external_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["oligo_name", "sequence_5to3", "orientation", "category", "target_name", "tm_c", "gc_percent"],
            )
            writer.writeheader()
            writer.writerows(external_rows)

        order_items = list(constant_order_items.values()) + external_rows + gene_specific_order_items
        with plate_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["plate", "well", "oligo_name", "target_name", "category", "orientation", "sequence_5to3", "stock_uM"]
            )
            for idx, item in enumerate(order_items):
                plate, well = _well_name(idx)
                writer.writerow(
                    [
                        plate,
                        well,
                        item["oligo_name"],
                        item["target_name"],
                        item["category"],
                        item["orientation"],
                        item["sequence_5to3"],
                        20,
                    ]
                )

        _write_fasta(order_fasta, ((item["oligo_name"], item["sequence_5to3"]) for item in order_items))
        instructions_txt.write_text(
            _render_assembly_instructions(
                target_designs=target_designs,
                constant_oligo_count=len(constant_order_items),
                external_forward_name="CONST_EXT_F",
                external_reverse_name="CONST_EXT_R",
            ),
            encoding="utf-8",
        )

        return {
            "oligos_csv": oligo_csv,
            "assembly_csv": assembly_csv,
            "constructs_fasta": constructs_fasta,
            "order_fasta": order_fasta,
            "plate_csv": plate_csv,
            "external_csv": external_csv,
            "instructions_txt": instructions_txt,
            "run_log": run_log,
        }
    finally:
        if owns_logger:
            for handler in managed_logger.handlers[:]:
                handler.close()
                managed_logger.removeHandler(handler)


def _default_length_scan_values(max_length: int) -> List[int]:
    anchors = {20, 30, 40, 50, 60, 80, 100, 120, 150, 200, max_length}
    values = sorted(value for value in anchors if MIN_OLIGO_LENGTH <= value <= max_length)
    if max_length not in values:
        values.append(max_length)
    return values


def run_length_optimize_gene_oligos(
    sequence_fasta: Path,
    output_dir: Path,
    input_type: str = "auto",
    max_length_cap: int = 120,
    candidate_lengths: Optional[Sequence[int]] = None,
    max_order_once_length: int = DEFAULT_MAX_ORDER_ONCE_LENGTH,
    target_host: str = DEFAULT_TARGET_HOST,
    tag_mode: str = "none",
    log_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    config: Optional[Dict[str, object]] = None,
) -> Dict[str, Path]:
    if max_length_cap < MIN_OLIGO_LENGTH:
        raise ValueError(f"max_length_cap must be at least {MIN_OLIGO_LENGTH}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_log = Path(log_path) if log_path else output_dir / "run.log"
    managed_logger, owns_logger = _setup_logger(run_log, logger)

    try:
        scan_values = list(candidate_lengths) if candidate_lengths else _default_length_scan_values(max_length_cap)
        scan_values = sorted({value for value in scan_values if value >= MIN_OLIGO_LENGTH and value <= max_length_cap})
        if not scan_values:
            raise ValueError("No candidate max-length values remain after filtering.")

        results: List[Dict[str, object]] = []
        for candidate in scan_values:
            with tempfile.TemporaryDirectory(prefix="uht_gene_opt_") as tmp_dir:
                tmp_out = Path(tmp_dir)
                try:
                    run_design_gene_oligos(
                        sequence_fasta=sequence_fasta,
                        output_dir=tmp_out,
                        input_type=input_type,
                        target_oligo_length=candidate,
                        max_order_once_length=max_order_once_length,
                        target_host=target_host,
                        tag_mode=tag_mode,
                        config=config,
                    )
                    assembly_rows = list(csv.DictReader((tmp_out / "assembly_report.csv").open()))
                    oligo_rows = list(csv.DictReader((tmp_out / "gene_oligos.csv").open()))
                    gene_specific_rows = [row for row in oligo_rows if row["oligo_role"] == "gene_specific"]
                    constant_rows = [row for row in oligo_rows if row["oligo_role"] == "constant"]
                    total_gene_specific = len(gene_specific_rows)
                    observed_lengths = [int(row["oligo_length"]) for row in gene_specific_rows]
                    max_observed = max(observed_lengths) if observed_lengths else 0
                    avg_length = (sum(observed_lengths) / len(observed_lengths)) if observed_lengths else 0.0
                    strategy_summary = "; ".join(
                        f"{row['target_name']}:{row['assembly_strategy']}" for row in assembly_rows
                    )
                    results.append(
                        {
                            "candidate_max_length_nt": candidate,
                            "status": "ok",
                            "gene_specific_primer_count": total_gene_specific,
                            "constant_primer_count": len(constant_rows),
                            "observed_max_gene_specific_length_nt": max_observed,
                            "observed_mean_gene_specific_length_nt": f"{avg_length:.2f}",
                            "strategy_summary": strategy_summary,
                            "error": "",
                        }
                    )
                except Exception as exc:
                    managed_logger.warning("Length candidate %s failed: %s", candidate, exc)
                    results.append(
                        {
                            "candidate_max_length_nt": candidate,
                            "status": "failed",
                            "gene_specific_primer_count": "",
                            "constant_primer_count": "",
                            "observed_max_gene_specific_length_nt": "",
                            "observed_mean_gene_specific_length_nt": "",
                            "strategy_summary": "",
                            "error": str(exc),
                        }
                    )

        successful = [row for row in results if row["status"] == "ok"]
        if not successful:
            raise ValueError("No scanned max-length values produced a valid design.")

        successful.sort(
            key=lambda row: (
                int(row["gene_specific_primer_count"]),
                int(row["observed_max_gene_specific_length_nt"]),
                float(row["observed_mean_gene_specific_length_nt"]),
                int(row["candidate_max_length_nt"]),
            )
        )
        best = successful[0]

        csv_path = output_dir / "length_optimization.csv"
        summary_path = output_dir / "length_optimization_summary.txt"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "candidate_max_length_nt",
                    "status",
                    "gene_specific_primer_count",
                    "constant_primer_count",
                    "observed_max_gene_specific_length_nt",
                    "observed_mean_gene_specific_length_nt",
                    "strategy_summary",
                    "error",
                ],
            )
            writer.writeheader()
            writer.writerows(results)

        summary_lines = [
            "Length optimization summary",
            "",
            f"Scanned max gene-specific oligo lengths: {', '.join(str(value) for value in scan_values)}",
            f"Recommended max length: {best['candidate_max_length_nt']} nt",
            f"Gene-specific primer count: {best['gene_specific_primer_count']}",
            f"Observed max gene-specific oligo length: {best['observed_max_gene_specific_length_nt']} nt",
            f"Observed mean gene-specific oligo length: {best['observed_mean_gene_specific_length_nt']} nt",
            f"Assembly strategies: {best['strategy_summary']}",
            "",
            "Ranking rule: minimize gene-specific primer count first, then prefer shorter resulting oligos.",
        ]
        summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

        return {
            "optimization_csv": csv_path,
            "optimization_summary": summary_path,
            "run_log": run_log,
        }
    finally:
        if owns_logger:
            for handler in managed_logger.handlers[:]:
                handler.close()
                managed_logger.removeHandler(handler)

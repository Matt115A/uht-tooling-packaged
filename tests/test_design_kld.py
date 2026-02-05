import csv
from pathlib import Path

import pytest
from Bio.Seq import Seq

from uht_tooling.workflows import design_kld


def write_fasta(path: Path, name: str, seq: str) -> None:
    path.write_text(f">{name}\n{seq}\n")


def test_degenerate_helpers():
    assert design_kld.is_valid_degenerate_codon("NNK")
    assert not design_kld.is_valid_degenerate_codon("AXT")
    assert not design_kld.is_valid_degenerate_codon("NN")

    assert design_kld.contains_degenerate_bases("ATGN")
    assert not design_kld.contains_degenerate_bases("ATGC")

    expanded = set(design_kld.expand_degenerate_sequence("RY"))
    assert expanded == {"AC", "AT", "GC", "GT"}


def test_pick_mutant_codon_min_diff():
    assert design_kld.pick_mutant_codon("GAA", "K") == "AAA"


def test_design_forward_reverse_primer(monkeypatch):
    monkeypatch.setattr(design_kld, "calc_tm_binding_region", lambda seq: 55.0)

    full_seq = "ATGCGT" * 30
    mutation_nt_pos = 30
    new_codon = "AAA"

    fwd_seq, fwd_binding, fwd_tm, fwd_gc, fwd_len = design_kld.design_forward_primer(
        full_seq, mutation_nt_pos, new_codon
    )
    assert fwd_seq.startswith(new_codon)
    assert len(fwd_binding) == design_kld.MIN_LENGTH
    assert fwd_len == len(new_codon) + design_kld.MIN_LENGTH
    assert fwd_tm == 55.0
    assert 40.0 <= fwd_gc <= 60.0

    rev_seq, rev_tm, rev_gc, rev_len = design_kld.design_reverse_primer(
        full_seq, mutation_nt_pos
    )
    upstream_seq = full_seq[mutation_nt_pos - design_kld.MIN_LENGTH : mutation_nt_pos]
    expected_rev = str(Seq(upstream_seq).reverse_complement())
    assert rev_seq == expected_rev
    assert rev_len == design_kld.MIN_LENGTH
    assert rev_tm == 55.0
    assert 40.0 <= rev_gc <= 60.0


def test_balance_primer_tms_trims_forward(monkeypatch):
    monkeypatch.setattr(design_kld, "calc_tm_binding_region", lambda seq: len(seq) + 40.0)

    full_seq = "ATGCGT" * 20
    mutation_nt_pos = 30
    new_codon = "AAA"

    fwd_binding = "ATGCGT" * 3 + "AT"  # length 20
    fwd_tm = design_kld.calc_tm_binding_region(fwd_binding)
    fwd_seq = new_codon + fwd_binding
    fwd_gc = design_kld.calc_gc_content(fwd_seq)
    fwd_result = (fwd_seq, fwd_binding, fwd_tm, fwd_gc, len(fwd_seq))

    rev_seq = "T" * 18
    rev_tm = 50.0
    rev_gc = design_kld.calc_gc_content(rev_seq)
    rev_result = (rev_seq, rev_tm, rev_gc, len(rev_seq))

    new_fwd, new_rev = design_kld.balance_primer_tms(
        fwd_result, rev_result, full_seq, mutation_nt_pos, new_codon
    )

    new_fwd_seq, new_fwd_binding, new_fwd_tm, _, _ = new_fwd
    assert len(new_fwd_binding) == 15
    assert new_fwd_tm == 55.0
    assert new_fwd_seq.startswith(new_codon)
    assert new_rev == rev_result


def test_run_design_kld_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(design_kld, "calc_tm_binding_region", lambda seq: 55.0)
    monkeypatch.setattr(design_kld, "calc_gc_content", lambda seq: 50.0)

    prefix = "AGCTTC" * 10
    gene = "ATGCGT" * 10
    suffix = "TCGATC" * 10
    context = prefix + gene + suffix

    gene_fasta = tmp_path / "gene.fasta"
    context_fasta = tmp_path / "context.fasta"
    mutations_csv = tmp_path / "mutations.csv"

    write_fasta(gene_fasta, "gene", gene)
    write_fasta(context_fasta, "context", context)

    with mutations_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mutations"])
        writer.writerow(["M1A"])
        writer.writerow(["R2:NNK"])
        writer.writerow(["M3MAA"])
        writer.writerow(["R4Del"])

    output_dir = tmp_path / "out"
    results_path = design_kld.run_design_kld(
        gene_fasta=gene_fasta,
        context_fasta=context_fasta,
        mutations_csv=mutations_csv,
        output_dir=output_dir,
    )

    assert results_path.exists()

    rows = list(csv.reader(results_path.open()))
    assert rows[0] == ["Primer Name", "Sequence", "Tm (binding)", "GC%", "Length", "Notes"]
    assert len(rows) == 1 + 2 * 4

    names = [row[0] for row in rows[1:]]
    assert "M1A_F" in names
    assert "M1A_R" in names
    assert "R2:NNK_F" in names
    assert "R4Del_F" in names

    notes = {row[0]: row[5] for row in rows[1:]}
    assert notes["R4Del_F"] == "Deletion"
    assert notes["M3MAA_F"] == "Insertion: +2 AA"
    assert notes["R2:NNK_F"].startswith("Library:")

    lengths = {row[0]: int(row[4]) for row in rows[1:]}
    assert lengths["M1A_F"] == 21
    assert lengths["R2:NNK_F"] == 21
    assert lengths["M3MAA_F"] == 24
    assert lengths["R4Del_F"] == 18


def test_run_design_kld_invalid_degenerate_codon(tmp_path, monkeypatch):
    monkeypatch.setattr(design_kld, "calc_tm_binding_region", lambda seq: 55.0)
    monkeypatch.setattr(design_kld, "calc_gc_content", lambda seq: 50.0)

    gene = "ATGCGT" * 10
    context = "AGCTTC" * 10 + gene + "TCGATC" * 10

    gene_fasta = tmp_path / "gene.fasta"
    context_fasta = tmp_path / "context.fasta"
    mutations_csv = tmp_path / "mutations.csv"

    write_fasta(gene_fasta, "gene", gene)
    write_fasta(context_fasta, "context", context)

    with mutations_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mutations"])
        writer.writerow(["R2:XXZ"])

    output_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Invalid degenerate codon"):
        design_kld.run_design_kld(
            gene_fasta=gene_fasta,
            context_fasta=context_fasta,
            mutations_csv=mutations_csv,
            output_dir=output_dir,
        )

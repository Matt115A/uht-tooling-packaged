import csv
import sys
import types
from pathlib import Path

from Bio.Seq import Seq
from typer.testing import CliRunner

from uht_tooling.workflows import design_synthetic_gene_pool


def install_fake_dnachisel(monkeypatch):
    module = types.ModuleType("dnachisel")

    class EnforceTranslation:
        def __init__(self, *args, **kwargs):
            pass

    class CodonOptimize:
        def __init__(self, species, location=None):
            self.species = species
            self.location = location

    class DnaOptimizationProblem:
        def __init__(self, sequence, constraints=None, objectives=None):
            self.sequence = sequence
            self.objectives = objectives or []

        def resolve_constraints(self):
            return None

        def optimize(self):
            if self.objectives and getattr(self.objectives[0], "species", "") == "s_cerevisiae":
                self.sequence = self.sequence.replace("GCT", "GCC")
            return None

    module.EnforceTranslation = EnforceTranslation
    module.CodonOptimize = CodonOptimize
    module.DnaOptimizationProblem = DnaOptimizationProblem
    monkeypatch.setitem(sys.modules, "dnachisel", module)


def write_fasta_records(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def base_config() -> dict:
    return {
        "gene_oligos": {
            "constant_5prime_dna": "TAATACGACTCACTATAGGGAGA",
            "constant_3prime_dna": "GCGTTTTTTTTTGGGGAAAA",
            "stop_codon": "TAA",
            "tags": {
                "n_his": {"dna": "CATCATCATCATCATCAT"},
                "c_his": {"dna": "CATCATCATCATCATCAT"},
                "n_his_flag": {"dna": "CATCATCATCATCATCATGACTACAAAGACGATGACGACAAG"},
                "c_his_flag": {"dna": "CATCATCATCATCATCATGACTACAAAGACGATGACGACAAG"},
            },
            "codon_table": {
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
            },
        }
    }


def test_run_design_synthetic_gene_pool_outputs(tmp_path, monkeypatch):
    install_fake_dnachisel(monkeypatch)

    input_fasta = tmp_path / "targets.fasta"
    write_fasta_records(input_fasta, [("gene_a", "ATGGCTGCTTAA"), ("gene_b", "AXK")])
    output_dir = tmp_path / "out"

    outputs = design_synthetic_gene_pool.run_design_synthetic_gene_pool(
        sequence_fasta=input_fasta,
        output_dir=output_dir,
        input_type="auto",
        tag_mode="n_his",
        config=base_config(),
    )

    for path in outputs.values():
        assert path.exists()

    pool_rows = list(csv.DictReader((output_dir / "synthetic_gene_pool.csv").open()))
    assert len(pool_rows) == 2
    assert all(row["pool_oligo_name"].endswith("_POOL") for row in pool_rows)
    gene_b_row = next(row for row in pool_rows if row["target_name"] == "gene_b")
    assert "NNN" in gene_b_row["sequence_5to3"]
    assert "Translated protein residue X to NNN" in gene_b_row["warnings"]
    primer_rows = list(csv.DictReader((output_dir / "synthetic_gene_pool_primers.csv").open()))
    assert [row["primer_name"] for row in primer_rows] == ["POOL_CONST_F", "POOL_CONST_R"]
    ordering_lines = (output_dir / "synthetic_gene_pool_ordering.tsv").read_text().splitlines()
    assert any("common_primer" in line for line in ordering_lines)


def test_design_synthetic_gene_pool_cli_smoke(tmp_path, monkeypatch):
    install_fake_dnachisel(monkeypatch)
    fake_align_apps = types.ModuleType("Bio.Align.Applications")
    fake_align_apps.MafftCommandline = object
    monkeypatch.setitem(sys.modules, "Bio.Align.Applications", fake_align_apps)
    from uht_tooling.cli import app

    input_fasta = tmp_path / "targets.fasta"
    write_fasta_records(input_fasta, [("gene_a", "AK"), ("gene_b", "MAK")])

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "design-synthetic-gene-pool",
            "--sequence-fasta",
            str(input_fasta),
            "--output-dir",
            str(tmp_path / "out"),
            "--input-type",
            "protein",
            "--tag-mode",
            "n_his",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "Synthetic gene pool outputs written" in result.stdout


def test_design_synthetic_gene_pool_enforces_common_primer_length_limit(tmp_path, monkeypatch):
    install_fake_dnachisel(monkeypatch)

    cfg = base_config()
    cfg["gene_oligos"]["constant_5prime_dna"] = "A" * 97
    input_fasta = tmp_path / "targets.fasta"
    write_fasta_records(input_fasta, [("gene_a", "AK")])

    try:
        design_synthetic_gene_pool.run_design_synthetic_gene_pool(
            sequence_fasta=input_fasta,
            output_dir=tmp_path / "out",
            input_type="protein",
            tag_mode="none",
            config=cfg,
        )
    except ValueError as exc:
        assert "99 bp ordering limit" in str(exc)
    else:
        raise AssertionError("Expected a primer-length limit error.")


def test_design_synthetic_gene_pool_pullout_primers_are_deterministic(tmp_path, monkeypatch):
    install_fake_dnachisel(monkeypatch)

    input_fasta_a = tmp_path / "targets_a.fasta"
    input_fasta_b = tmp_path / "targets_b.fasta"
    write_fasta_records(input_fasta_a, [("gene_a", "AK"), ("gene_b", "MAK")])
    write_fasta_records(input_fasta_b, [("gene_x", "MK"), ("gene_y", "AAK")])

    output_dir_a = tmp_path / "out_a"
    output_dir_b = tmp_path / "out_b"

    design_synthetic_gene_pool.run_design_synthetic_gene_pool(
        sequence_fasta=input_fasta_a,
        output_dir=output_dir_a,
        input_type="protein",
        tag_mode="none",
        include_pullout_primers=True,
        config=base_config(),
    )
    design_synthetic_gene_pool.run_design_synthetic_gene_pool(
        sequence_fasta=input_fasta_b,
        output_dir=output_dir_b,
        input_type="protein",
        tag_mode="none",
        include_pullout_primers=True,
        config=base_config(),
    )

    pool_rows_a = list(csv.DictReader((output_dir_a / "synthetic_gene_pool.csv").open()))
    pool_rows_b = list(csv.DictReader((output_dir_b / "synthetic_gene_pool.csv").open()))
    assert [row["pullout_index"] for row in pool_rows_a] == ["AACAACAACC", "ACATAACACC"]
    assert [row["pullout_index"] for row in pool_rows_b] == ["AACAACAACC", "ACATAACACC"]
    assert pool_rows_a[0]["sequence_5to3"].endswith("TAAAACAACAACCGCGTTTTTTTTTGGGGAAAA")
    assert pool_rows_a[1]["sequence_5to3"].endswith("TAAACATAACACCGCGTTTTTTTTTGGGGAAAA")

    primer_rows = list(csv.DictReader((output_dir_a / "synthetic_gene_pool_primers.csv").open()))
    common_names = [row["primer_name"] for row in primer_rows[:2]]
    assert common_names == ["POOL_CONST_F", "POOL_CONST_R"]
    pullout_rows = [row for row in primer_rows if row["primer_role"] == "gene_specific_pullout"]
    assert [row["primer_name"] for row in pullout_rows] == ["gene_a_PULLOUT_R", "gene_b_PULLOUT_R"]
    assert [row["pullout_index"] for row in pullout_rows] == ["AACAACAACC", "ACATAACACC"]
    assert all(row["adds_sequence"] == "" for row in pullout_rows)
    assert pullout_rows[0]["anneal_sequence"] == str(Seq("AACAACAACCGCGTTTTTTTTTGGGGAAAA").reverse_complement())
    assert pullout_rows[1]["anneal_sequence"] == str(Seq("ACATAACACCGCGTTTTTTTTTGGGGAAAA").reverse_complement())

    ordering_lines = (output_dir_a / "synthetic_gene_pool_ordering.tsv").read_text().splitlines()
    assert any("pullout_primer" in line for line in ordering_lines)


def test_pullout_index_table_metadata_is_packaged():
    table = design_synthetic_gene_pool._load_pullout_index_table()

    assert table["table_name"] == "cfps_pullout_indexes"
    assert table["entry_count"] >= 100
    assert table["selection_constraints"]["length_nt"] == 10
    assert table["selection_constraints"]["min_pairwise_hamming_distance"] == 5
    assert table["table_version"] == "2026-06-15-v2"

    first = table["entries"][0]
    assert first["sequence"] == "AACAACAACC"
    assert 40.0 <= float(first["gc_percent"]) <= 60.0
    assert int(first["min_hamming_distance_to_table"]) >= 5

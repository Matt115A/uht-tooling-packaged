import csv
import sys
import types
from pathlib import Path

from typer.testing import CliRunner

from uht_tooling.workflows import design_gene_oligos


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
            self.constraints = constraints or []
            self.objectives = objectives or []

        def resolve_constraints(self):
            return None

        def optimize(self):
            if self.objectives:
                species = getattr(self.objectives[0], "species", "")
                if species == "s_cerevisiae":
                    self.sequence = self.sequence.replace("GCT", "GCC")
                elif species == "h_sapiens":
                    self.sequence = self.sequence.replace("AAA", "AAG")
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


def test_detect_input_type_and_reverse_translate():
    assert design_gene_oligos._detect_input_type("ATGGCC") == "dna"
    assert design_gene_oligos._detect_input_type("MKT") == "protein"
    assert design_gene_oligos.reverse_translate_protein(
        "MAK", {"M": "ATG", "A": "GCT", "K": "AAA"}
    ) == "ATGGCTAAA"


def test_run_design_gene_oligos_multi_target_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)

    input_fasta = tmp_path / "targets.fasta"
    write_fasta_records(
        input_fasta,
        [
            ("gene_a", "GCTGCTGCTGCTGCTGCTTAA"),
            ("gene_b", "ATGGCTGCTGCTGCTGCTTAA"),
        ],
    )

    output_dir = tmp_path / "out"
    outputs = design_gene_oligos.run_design_gene_oligos(
        sequence_fasta=input_fasta,
        output_dir=output_dir,
        input_type="dna",
        target_oligo_length=40,
        tag_mode="none",
        config=base_config(),
    )

    for path in outputs.values():
        assert path.exists()

    oligo_rows = list(csv.DictReader((output_dir / "gene_oligos.csv").open()))
    assert {row["target_name"] for row in oligo_rows} == {"ALL", "gene_a", "gene_b"}
    assert any(row["oligo_role"] == "constant" for row in oligo_rows)
    assert any(row["oligo_role"] == "gene_specific" for row in oligo_rows)

    assembly_rows = list(csv.DictReader((output_dir / "assembly_report.csv").open()))
    assembly_by_name = {row["target_name"]: row for row in assembly_rows}
    assert assembly_by_name["gene_a"]["added_start_codon"] == "yes"
    assert assembly_by_name["gene_b"]["added_start_codon"] == "no"
    assert all(row["added_stop_codon"] == "no" for row in assembly_rows)

    plate_rows = list(csv.DictReader((output_dir / "ordering_plate.csv").open()))
    assert any(row["category"] in {"constant", "external"} and row["target_name"] == "ALL" for row in plate_rows)
    assert any(row["category"] == "gene_specific" and row["target_name"] == "gene_a" for row in plate_rows)
    assert any(row["category"] == "external" for row in plate_rows)


def test_run_design_gene_oligos_n_terminal_tag_gets_start_codon(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    install_fake_dnachisel(monkeypatch)

    protein_fasta = tmp_path / "protein.fasta"
    write_fasta_records(protein_fasta, [("protein", "AK*")])

    output_dir = tmp_path / "out"
    design_gene_oligos.run_design_gene_oligos(
        sequence_fasta=protein_fasta,
        output_dir=output_dir,
        input_type="protein",
        target_oligo_length=40,
        tag_mode="n_his",
        config=base_config(),
    )

    fasta_lines = (output_dir / "final_constructs.fasta").read_text().splitlines()
    construct_seq = "".join(line for line in fasta_lines[1:] if not line.startswith(">"))
    cfg = base_config()["gene_oligos"]
    assert construct_seq.startswith(cfg["constant_5prime_dna"] + "ATG" + cfg["tags"]["n_his"]["dna"])
    assert "GCTAAA" in construct_seq


def test_run_design_gene_oligos_uses_builtin_defaults_without_config(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    install_fake_dnachisel(monkeypatch)

    protein_fasta = tmp_path / "protein.fasta"
    write_fasta_records(protein_fasta, [("protein", "MAK")])

    output_dir = tmp_path / "out"
    design_gene_oligos.run_design_gene_oligos(
        sequence_fasta=protein_fasta,
        output_dir=output_dir,
        input_type="protein",
        target_oligo_length=40,
        tag_mode="n_his",
    )

    construct_seq = "".join(
        line
        for line in (output_dir / "final_constructs.fasta").read_text().splitlines()
        if not line.startswith(">")
    )
    assert construct_seq.startswith(design_gene_oligos.DEFAULT_CONSTANT_5PRIME_DNA)
    assert construct_seq.endswith(design_gene_oligos.DEFAULT_CONSTANT_3PRIME_DNA)


def test_tile_sequence_warns_when_overlap_outside_window(monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 70.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)

    tiles = design_gene_oligos.tile_sequence_with_overlaps("A" * 90, 40)

    assert any(tile.warning for tile in tiles[1:])
    assert all(20 <= tile.length <= 50 for tile in tiles)


def test_large_requested_oligo_length_is_treated_as_maximum(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)

    input_fasta = tmp_path / "target.fasta"
    write_fasta_records(input_fasta, [("gene_a", "ATG" + "GCT" * 25 + "TAA")])

    output_dir = tmp_path / "out"
    design_gene_oligos.run_design_gene_oligos(
        sequence_fasta=input_fasta,
        output_dir=output_dir,
        input_type="dna",
        target_oligo_length=200,
        tag_mode="none",
        config=base_config(),
    )

    rows = [row for row in csv.DictReader((output_dir / "gene_oligos.csv").open()) if row["target_name"] == "gene_a"]
    assert rows
    assert max(int(row["oligo_length"]) for row in rows) <= 200


def test_design_gene_oligos_cli_smoke(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    fake_align_apps = types.ModuleType("Bio.Align.Applications")
    fake_align_apps.MafftCommandline = object
    monkeypatch.setitem(sys.modules, "Bio.Align.Applications", fake_align_apps)
    install_fake_dnachisel(monkeypatch)
    from uht_tooling.cli import app

    input_fasta = tmp_path / "protein.fasta"
    write_fasta_records(input_fasta, [("protein_a", "AK"), ("protein_b", "MAK")])

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "design-gene-oligos",
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
    assert "Gene oligo outputs written" in result.stdout


def test_length_optimize_gene_oligos_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    install_fake_dnachisel(monkeypatch)

    input_fasta = tmp_path / "target.fasta"
    write_fasta_records(input_fasta, [("gene_a", "ATG" + "GCT" * 30 + "TAA")])

    output_dir = tmp_path / "opt_out"
    outputs = design_gene_oligos.run_length_optimize_gene_oligos(
        sequence_fasta=input_fasta,
        output_dir=output_dir,
        input_type="dna",
        max_length_cap=120,
        candidate_lengths=[30, 45, 80],
        tag_mode="none",
        config=base_config(),
    )

    for path in outputs.values():
        assert path.exists()

    rows = list(csv.DictReader((output_dir / "length_optimization.csv").open()))
    assert [int(row["candidate_max_length_nt"]) for row in rows] == [30, 45, 80]
    assert all(row["status"] == "ok" for row in rows)
    summary = (output_dir / "length_optimization_summary.txt").read_text()
    assert "Recommended max length:" in summary


def test_length_optimize_gene_oligos_cli_smoke(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    fake_align_apps = types.ModuleType("Bio.Align.Applications")
    fake_align_apps.MafftCommandline = object
    monkeypatch.setitem(sys.modules, "Bio.Align.Applications", fake_align_apps)
    install_fake_dnachisel(monkeypatch)
    from uht_tooling.cli import app

    input_fasta = tmp_path / "protein.fasta"
    write_fasta_records(input_fasta, [("protein_a", "AK"), ("protein_b", "MAK")])

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "length-optimize-gene-oligos",
            "--sequence-fasta",
            str(input_fasta),
            "--output-dir",
            str(tmp_path / "opt_out"),
            "--input-type",
            "protein",
            "--scan-length",
            "30",
            "--scan-length",
            "45",
            "--tag-mode",
            "n_his",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "Length optimization outputs written" in result.stdout


def test_protein_input_uses_target_host_for_codon_optimization(tmp_path, monkeypatch):
    monkeypatch.setattr(design_gene_oligos, "calc_tm", lambda seq: 60.0)
    monkeypatch.setattr(design_gene_oligos, "calc_gc_content", lambda seq: 50.0)
    install_fake_dnachisel(monkeypatch)

    input_fasta = tmp_path / "protein.fasta"
    write_fasta_records(input_fasta, [("protein", "AK")])

    output_dir = tmp_path / "out"
    design_gene_oligos.run_design_gene_oligos(
        sequence_fasta=input_fasta,
        output_dir=output_dir,
        input_type="protein",
        target_oligo_length=40,
        target_host="s_cerevisiae",
        tag_mode="none",
        config=base_config(),
    )

    construct_seq = "".join(
        line for line in (output_dir / "final_constructs.fasta").read_text().splitlines() if not line.startswith(">")
    )
    assert "GCCAAA" in construct_seq

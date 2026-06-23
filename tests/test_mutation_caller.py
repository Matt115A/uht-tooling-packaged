from types import SimpleNamespace

from uht_tooling.workflows import mutation_caller


def test_align_to_reference_preserves_original_fastq_headers(monkeypatch):
    gene_reads = {
        "7feee2ef-7cb2-460e-8e09-e05f65a1f835 runid=abc read=1": ("ACGT", "IIII"),
    }

    class FakePopen:
        def __init__(self, *args, **kwargs):
            self.returncode = 0

        def communicate(self):
            return ">REF\nACGT\n>READ_0\nACGT\n", ""

    monkeypatch.setattr(mutation_caller, "MafftCommandline", lambda input: "mafft fake")
    monkeypatch.setattr(mutation_caller.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(
        mutation_caller.AlignIO,
        "read",
        lambda path, fmt: [
            SimpleNamespace(id="REF", seq="ACGT"),
            SimpleNamespace(id="READ_0", seq="ACGT"),
        ],
    )

    aligned_ref, aligned_reads, aligned_quals = mutation_caller.align_to_reference(gene_reads, "ACGT")

    original_id = next(iter(gene_reads))
    assert aligned_ref == "ACGT"
    assert aligned_reads == {original_id: "ACGT"}
    assert aligned_quals == {original_id: "IIII"}

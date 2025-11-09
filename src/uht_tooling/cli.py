from pathlib import Path
from typing import Optional

import typer

from uht_tooling.workflows.design_gibson import run_design_gibson
from uht_tooling.workflows.design_slim import run_design_slim
from uht_tooling.workflows.mutation_caller import (
    expand_fastq_inputs as expand_fastq_inputs_mutation,
    run_mutation_caller,
)
from uht_tooling.workflows.nextera_designer import run_nextera_primer_design
from uht_tooling.workflows.profile_inserts import (
    expand_fastq_inputs as expand_fastq_inputs_profile,
    run_profile_inserts,
)
from uht_tooling.workflows.umi_hunter import (
    expand_fastq_inputs as expand_fastq_inputs_umi,
    run_umi_hunter,
)
from uht_tooling.workflows.gui import launch_gui

app = typer.Typer(help="Command-line interface for the uht-tooling package.")


@app.command("design-slim", help="Design SLIM primers from user-specified FASTA/CSV inputs.")
def design_slim_command(
    gene_fasta: Path = typer.Option(..., exists=True, readable=True, help="Path to the gene FASTA file."),
    context_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="Path to the context FASTA file containing the plasmid or genomic sequence.",
    ),
    mutations_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV file containing a 'mutations' column with the desired edits.",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory where results will be written.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file for this run.",
    ),
):
    """Design SLIM primers from user-provided inputs."""
    run_design_slim(
        gene_fasta=gene_fasta,
        context_fasta=context_fasta,
        mutations_csv=mutations_csv,
        output_dir=output_dir,
        log_path=log_path,
    )
    typer.echo(f"SLIM primers written to {output_dir / 'SLIM_primers.csv'}")


@app.command("nextera-primers", help="Generate Nextera XT primers from binding region CSV input.")
def nextera_primers_command(
    binding_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV file with a 'binding_region' column; first row is i7, second row is i5.",
    ),
    output_csv: Path = typer.Option(
        ...,
        dir_okay=False,
        writable=True,
        help="Path to write the generated primer CSV.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
    config: Optional[Path] = typer.Option(
        None,
        exists=True,
        readable=True,
        help="Optional YAML file providing overrides for indexes/prefixes/suffixes.",
    ),
):
    """Generate Nextera XT primers from user-supplied binding regions."""
    result_path = run_nextera_primer_design(
        binding_csv=binding_csv,
        output_csv=output_csv,
        log_path=log_path,
        config_path=config,
    )
    typer.echo(f"Nextera primers written to {result_path}")


@app.command("design-gibson", help="Design Gibson assembly primers and assembly plans.")
def design_gibson_command(
    gene_fasta: Path = typer.Option(..., exists=True, readable=True, help="Path to the gene FASTA file."),
    context_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="Path to the circular context FASTA file.",
    ),
    mutations_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV file with a 'mutations' column (use '+' to link sub-mutations).",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory where primer and assembly plan CSVs will be written.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path for a dedicated log file.",
    ),
):
    """Design Gibson assembly primers for user-defined mutations."""
    outputs = run_design_gibson(
        gene_fasta=gene_fasta,
        context_fasta=context_fasta,
        mutations_csv=mutations_csv,
        output_dir=output_dir,
        log_path=log_path,
    )
    typer.echo("Gibson outputs written:")
    for name, path in outputs.items():
        typer.echo(f"  {name}: {path}")


@app.command(
    "mutation-caller",
    help="Identify amino-acid substitutions from long-read data without UMIs.",
)
def mutation_caller_command(
    template_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="FASTA file containing the mutation caller template sequence.",
    ),
    flanks_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV file describing gene flanks and min/max lengths.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        help="One or more FASTQ(.gz) paths or glob patterns (provide multiple --fastq options as needed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory where per-sample outputs will be written.",
    ),
    threshold: int = typer.Option(
        10,
        min=1,
        help="Minimum AA substitution count to include in the frequent-substitution report.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
):
    """Identify and summarise amino-acid substitutions."""
    fastq_files = expand_fastq_inputs_mutation(fastq)
    results = run_mutation_caller(
        template_fasta=template_fasta,
        flanks_csv=flanks_csv,
        fastq_files=fastq_files,
        output_dir=output_dir,
        threshold=threshold,
        log_path=log_path,
    )
    if not results:
        typer.echo("No outputs were generated. Check inputs and threshold settings.")
    else:
        typer.echo("Mutation caller outputs:")
        for entry in results:
            typer.echo(f"  Sample {entry['sample']}: {entry['directory']}")


@app.command("umi-hunter", help="Cluster UMIs and produce consensus genes from long-read data.")
def umi_hunter_command(
    template_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="Template FASTA file for consensus generation.",
    ),
    config_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV describing UMI/gene flanks and length bounds.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory where UMI hunter outputs will be stored.",
    ),
    umi_identity_threshold: float = typer.Option(
        0.9,
        min=0.0,
        max=1.0,
        help="UMI clustering identity threshold (default: 0.9).",
    ),
    consensus_mutation_threshold: float = typer.Option(
        0.7,
        min=0.0,
        max=1.0,
        help="Mutation threshold for consensus calling (default: 0.7).",
    ),
    min_cluster_size: int = typer.Option(
        1,
        min=1,
        help="Minimum number of reads required in a UMI cluster before a consensus is generated.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
):
    """Cluster UMIs and generate consensus sequences from long-read FASTQ data."""
    fastq_files = expand_fastq_inputs_umi(fastq)
    results = run_umi_hunter(
        template_fasta=template_fasta,
        config_csv=config_csv,
        fastq_files=fastq_files,
        output_dir=output_dir,
        umi_identity_threshold=umi_identity_threshold,
        consensus_mutation_threshold=consensus_mutation_threshold,
        min_cluster_size=min_cluster_size,
        log_path=log_path,
    )
    if not results:
        typer.echo("No UMI hunter outputs generated.")
    else:
        typer.echo("UMI hunter outputs:")
        for entry in results:
            total_clusters = entry.get("clusters_total", entry.get("clusters", 0))
            typer.echo(
                f"  Sample {entry['sample']}: "
                f"{entry.get('clusters', 0)} consensus clusters "
                f"(from {total_clusters} total) → {entry['directory']}"
            )


@app.command("ep-library-profile", help="Profile mutation rates for ep-library sequencing data.")
def ep_library_profile_command(
    region_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="FASTA file describing the region of interest.",
    ),
    plasmid_fasta: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="FASTA file with the full plasmid sequence.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory for per-sample outputs.",
    ),
    work_dir: Optional[Path] = typer.Option(
        None,
        dir_okay=True,
        writable=True,
        help="Optional scratch directory for intermediate files (defaults to output/tmp).",
    ),
):
    """Quantify mutation rates for ep-library sequencing experiments."""
    fastq_files = expand_fastq_inputs_ep(fastq)
    results = run_ep_library_profile(
        fastq_paths=fastq_files,
        region_fasta=region_fasta,
        plasmid_fasta=plasmid_fasta,
        output_dir=output_dir,
        work_dir=work_dir,
    )
    samples = results.get("samples", [])
    if not samples:
        typer.echo("No ep-library profile outputs generated.")
    else:
        typer.echo(f"Master summary written to {results['master_summary']}")
        for sample in samples:
            typer.echo(f"  Sample {sample['sample']}: {sample['results_dir']}")


@app.command("profile-inserts", help="Extract and profile inserts using probe pairs.")
def profile_inserts_command(
    probes_csv: Path = typer.Option(
        ...,
        exists=True,
        readable=True,
        help="CSV file containing upstream/downstream probes.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        dir_okay=True,
        writable=True,
        help="Directory for per-sample outputs.",
    ),
    min_ratio: int = typer.Option(
        80,
        min=0,
        max=100,
        help="Minimum fuzzy match ratio for probe detection (default: 80).",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
):
    """Profile inserts in FASTQ reads using probe pairs and produce QC outputs."""
    fastq_files = expand_fastq_inputs_profile(fastq)
    results = run_profile_inserts(
        probes_csv=probes_csv,
        fastq_files=fastq_files,
        output_dir=output_dir,
        min_ratio=min_ratio,
        log_path=log_path,
    )
    if not results:
        typer.echo("No profile inserts outputs generated.")
    else:
        typer.echo("Profile inserts outputs:")
        for entry in results:
            typer.echo(f"  Sample {entry['sample']}: {entry['directory']}")


@app.command("gui", help="Launch the graphical interface.")
def gui_command(
    server_name: str = typer.Option(
        "127.0.0.1",
        "--server-name",
        "-n",
        help="Hostname or IP address to bind the GUI server.",
    ),
    server_port: Optional[int] = typer.Option(
        7860,
        "--server-port",
        "-p",
        help="Preferred port for the GUI (falls back automatically if unavailable).",
    ),
    share: bool = typer.Option(
        False,
        "--share",
        help="Enable Gradio's public sharing tunnel (requires network access).",
    ),
):
    """Launch the Gradio GUI."""
    try:
        launch_gui(server_name=server_name, server_port=server_port, share=share)
    except KeyboardInterrupt:
        typer.echo("GUI stopped by user.")
    except Exception as exc:
        typer.echo(f"Failed to start GUI: {exc}")
        raise typer.Exit(1)


def main():
    app()


if __name__ == "__main__":
    main()

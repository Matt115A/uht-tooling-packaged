from pathlib import Path
from typing import Optional

import typer

from uht_tooling.config import get_option, load_config
from uht_tooling.tools import ToolNotFoundError, validate_workflow_tools
from uht_tooling.workflows.design_gene_oligos import (
    DEFAULT_MAX_ORDER_ONCE_LENGTH,
    DEFAULT_TARGET_HOST,
    run_design_gene_oligos,
    run_length_optimize_gene_oligos,
)
from uht_tooling.workflows.design_synthetic_gene_pool import run_design_synthetic_gene_pool
from uht_tooling.workflows.design_gibson import run_design_gibson
from uht_tooling.workflows.design_kld import run_design_kld
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
from uht_tooling.workflows.mut_rate import (
    expand_fastq_inputs as expand_fastq_inputs_ep,
    run_ep_library_profile,
)
from uht_tooling.workflows.ssm_profiler import (
    expand_fastq_inputs as expand_fastq_inputs_ssm,
    parse_scheme_map,
    run_ssm_profiler,
)
from uht_tooling.workflows.enrich_monitor import (
    expand_fastq_inputs as expand_fastq_inputs_enrich,
    run_enrich_monitor,
)
from uht_tooling.web import launch_web_gui

app = typer.Typer(help="Command-line interface for the uht-tooling package.")


@app.callback()
def main_callback(
    ctx: typer.Context,
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-K",
        exists=True,
        readable=True,
        help="Path to YAML configuration file for default options.",
    ),
):
    """Global callback to load configuration file."""
    ctx.ensure_object(dict)
    ctx.obj["config"] = load_config(config)


@app.command("design-slim", help="Design SLIM primers from user-specified FASTA/CSV inputs.")
def design_slim_command(
    ctx: typer.Context,
    gene_fasta: Path = typer.Option(
        ..., "--gene-fasta", "-g", exists=True, readable=True, help="Path to the gene FASTA file."
    ),
    context_fasta: Path = typer.Option(
        ...,
        "--context-fasta",
        "-c",
        exists=True,
        readable=True,
        help="Path to the context FASTA file containing the plasmid or genomic sequence.",
    ),
    mutations_csv: Path = typer.Option(
        ...,
        "--mutations-csv",
        "-m",
        exists=True,
        readable=True,
        help="CSV file containing a 'mutations' column with the desired edits.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where results will be written.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
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


@app.command("design-kld", help="Design KLD (inverse PCR) primers from user-specified FASTA/CSV inputs.")
def design_kld_command(
    ctx: typer.Context,
    gene_fasta: Path = typer.Option(
        ..., "--gene-fasta", "-g", exists=True, readable=True, help="Path to the gene FASTA file."
    ),
    context_fasta: Path = typer.Option(
        ...,
        "--context-fasta",
        "-c",
        exists=True,
        readable=True,
        help="Path to the context FASTA file containing the plasmid or genomic sequence.",
    ),
    mutations_csv: Path = typer.Option(
        ...,
        "--mutations-csv",
        "-m",
        exists=True,
        readable=True,
        help="CSV file containing a 'mutations' column with the desired edits.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where results will be written.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file for this run.",
    ),
):
    """Design KLD (inverse PCR) primers from user-provided inputs."""
    result_path = run_design_kld(
        gene_fasta=gene_fasta,
        context_fasta=context_fasta,
        mutations_csv=mutations_csv,
        output_dir=output_dir,
        log_path=log_path,
    )
    typer.echo(f"KLD primers written to {result_path}")


@app.command("nextera-primers", help="Generate Nextera XT primers from binding region CSV input.")
def nextera_primers_command(
    ctx: typer.Context,
    binding_csv: Path = typer.Option(
        ...,
        "--binding-csv",
        "-b",
        exists=True,
        readable=True,
        help="CSV file with a 'binding_region' column; first row is i7, second row is i5.",
    ),
    output_csv: Path = typer.Option(
        ...,
        "--output-csv",
        "-o",
        dir_okay=False,
        writable=True,
        help="Path to write the generated primer CSV.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
    config: Optional[Path] = typer.Option(
        None,
        "--nextera-config",
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
    ctx: typer.Context,
    gene_fasta: Path = typer.Option(
        ..., "--gene-fasta", "-g", exists=True, readable=True, help="Path to the gene FASTA file."
    ),
    context_fasta: Path = typer.Option(
        ...,
        "--context-fasta",
        "-c",
        exists=True,
        readable=True,
        help="Path to the circular context FASTA file.",
    ),
    mutations_csv: Path = typer.Option(
        ...,
        "--mutations-csv",
        "-m",
        exists=True,
        readable=True,
        help="CSV file with a 'mutations' column (use '+' to link sub-mutations).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where primer and assembly plan CSVs will be written.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
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


@app.command("design-gene-oligos", help="Design overlap-extension PCR oligos for IVTT-ready gene constructs.")
def design_gene_oligos_command(
    ctx: typer.Context,
    sequence_fasta: Path = typer.Option(
        ..., "--sequence-fasta", "-s", exists=True, readable=True, help="Path to the DNA or protein FASTA file."
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where oligo design outputs will be written.",
    ),
    input_type: str = typer.Option(
        "auto",
        "--input-type",
        help="Interpret the input as auto, dna, or protein.",
    ),
    target_oligo_length: int = typer.Option(
        40,
        "--target-oligo-length",
        min=20,
        help="Maximum allowed gene-specific oligo length in nucleotides.",
    ),
    max_order_once_length: int = typer.Option(
        DEFAULT_MAX_ORDER_ONCE_LENGTH,
        "--max-orderonce-length",
        min=20,
        help="Maximum allowed length for reusable order-once oligos.",
    ),
    target_host: str = typer.Option(
        DEFAULT_TARGET_HOST,
        "--target-host",
        help="Target host used for automatic codon optimization of protein inputs.",
    ),
    tag_mode: str = typer.Option(
        "none",
        "--tag-mode",
        help="One of: none, n_his, c_his, n_his_flag, c_his_flag.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file for this run.",
    ),
):
    """Design overlap-extension PCR oligos for a single IVTT-ready construct."""
    config = ctx.obj.get("config", {}) if ctx.obj else {}
    target_oligo_length = get_option(
        config,
        "target_oligo_length",
        target_oligo_length,
        default=40,
        workflow="design_gene_oligos",
    )
    outputs = run_design_gene_oligos(
        sequence_fasta=sequence_fasta,
        output_dir=output_dir,
        input_type=input_type,
        target_oligo_length=int(target_oligo_length),
        max_order_once_length=int(max_order_once_length),
        target_host=target_host,
        tag_mode=tag_mode,
        log_path=log_path,
        config=config,
    )
    typer.echo("Gene oligo outputs written:")
    for name, path in outputs.items():
        typer.echo(f"  {name}: {path}")


@app.command("length-optimize-gene-oligos", help="Scan max-length caps and rank gene-oligo assembly strategies.")
def length_optimize_gene_oligos_command(
    ctx: typer.Context,
    sequence_fasta: Path = typer.Option(
        ..., "--sequence-fasta", "-s", exists=True, readable=True, help="Path to the DNA or protein FASTA file."
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where optimization outputs will be written.",
    ),
    input_type: str = typer.Option(
        "auto",
        "--input-type",
        help="Interpret the input as auto, dna, or protein.",
    ),
    max_length_cap: int = typer.Option(
        120,
        "--max-length-cap",
        min=20,
        help="Upper bound for scanned gene-specific oligo lengths.",
    ),
    scan_length: list[int] = typer.Option(
        [],
        "--scan-length",
        help="Explicit max-length value(s) to test. Repeat as needed; defaults to a built-in scan set.",
    ),
    max_order_once_length: int = typer.Option(
        DEFAULT_MAX_ORDER_ONCE_LENGTH,
        "--max-orderonce-length",
        min=20,
        help="Maximum allowed length for reusable order-once oligos.",
    ),
    target_host: str = typer.Option(
        DEFAULT_TARGET_HOST,
        "--target-host",
        help="Target host used for automatic codon optimization of protein inputs.",
    ),
    tag_mode: str = typer.Option(
        "none",
        "--tag-mode",
        help="One of: none, n_his, c_his, n_his_flag, c_his_flag.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file for this run.",
    ),
):
    """Scan candidate max-length caps and recommend a primer-economical strategy."""
    config = ctx.obj.get("config", {}) if ctx.obj else {}
    outputs = run_length_optimize_gene_oligos(
        sequence_fasta=sequence_fasta,
        output_dir=output_dir,
        input_type=input_type,
        max_length_cap=max_length_cap,
        candidate_lengths=scan_length or None,
        max_order_once_length=max_order_once_length,
        target_host=target_host,
        tag_mode=tag_mode,
        log_path=log_path,
        config=config,
    )
    typer.echo("Length optimization outputs written:")
    for name, path in outputs.items():
        typer.echo(f"  {name}: {path}")


@app.command("design-synthetic-gene-pool", help="Design pooled synthetic-gene ordering oligos and common lift-out primers.")
def design_synthetic_gene_pool_command(
    ctx: typer.Context,
    sequence_fasta: Path = typer.Option(
        ..., "--sequence-fasta", "-s", exists=True, readable=True, help="Path to the DNA or protein FASTA file."
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where synthetic-gene pool outputs will be written.",
    ),
    input_type: str = typer.Option(
        "auto",
        "--input-type",
        help="Interpret the input as auto, dna, or protein.",
    ),
    tag_mode: str = typer.Option(
        "none",
        "--tag-mode",
        help="One of: none, n_his, c_his, n_his_flag, c_his_flag.",
    ),
    include_pullout_primers: bool = typer.Option(
        False,
        "--include-pullout-primers/--no-include-pullout-primers",
        help="Add one gene-specific reverse pullout primer per design; each ordered gene carries a deterministic unique index just upstream of the terminator, and the pullout primer binds that built-in region after whole-pool amplification.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file for this run.",
    ),
):
    """Design a pooled synthetic-gene ordering set with reusable primers."""
    config = ctx.obj.get("config", {}) if ctx.obj else {}
    outputs = run_design_synthetic_gene_pool(
        sequence_fasta=sequence_fasta,
        output_dir=output_dir,
        input_type=input_type,
        tag_mode=tag_mode,
        include_pullout_primers=include_pullout_primers,
        log_path=log_path,
        config=config,
    )
    typer.echo("Synthetic gene pool outputs written:")
    for name, path in outputs.items():
        typer.echo(f"  {name}: {path}")


@app.command(
    "mutation-caller",
    help="Identify amino-acid substitutions from long-read data without UMIs.",
)
def mutation_caller_command(
    ctx: typer.Context,
    template_fasta: Path = typer.Option(
        ...,
        "--template-fasta",
        "-t",
        exists=True,
        readable=True,
        help="FASTA file containing the mutation caller template sequence.",
    ),
    flanks_csv: Path = typer.Option(
        ...,
        "--flanks-csv",
        "-f",
        exists=True,
        readable=True,
        help="CSV file describing gene flanks and min/max lengths.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help="One or more FASTQ(.gz) paths or glob patterns (provide multiple --fastq options as needed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where per-sample outputs will be written.",
    ),
    threshold: int = typer.Option(
        10,
        "--threshold",
        "-T",
        min=1,
        help="Minimum AA substitution count to include in the frequent-substitution report.",
    ),
    min_flank_ratio: int = typer.Option(
        80,
        "--min-flank-ratio",
        "-r",
        min=0,
        max=100,
        help="Minimum fuzzy match ratio (0-100) for flank detection during read extraction.",
    ),
    min_base_qual: int = typer.Option(
        0,
        "--min-base-qual",
        "-Q",
        min=0,
        max=93,
        help=(
            "Minimum per-base Phred quality score required at every position of a codon "
            "for a substitution call at that codon (0 disables filtering)."
        ),
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
):
    """Identify and summarise amino-acid substitutions."""
    # Validate required external tools
    try:
        validate_workflow_tools("mutation_caller")
    except ToolNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    # Apply config defaults
    config = ctx.obj.get("config", {}) if ctx.obj else {}
    threshold = get_option(config, "threshold", threshold, default=10, workflow="mutation_caller")
    min_flank_ratio = get_option(
        config, "min_flank_ratio", min_flank_ratio, default=80, workflow="mutation_caller"
    )
    min_base_qual = get_option(
        config, "min_base_qual", min_base_qual, default=0, workflow="mutation_caller"
    )

    fastq_files = expand_fastq_inputs_mutation(fastq)
    results = run_mutation_caller(
        template_fasta=template_fasta,
        flanks_csv=flanks_csv,
        fastq_files=fastq_files,
        output_dir=output_dir,
        threshold=threshold,
        min_flank_ratio=min_flank_ratio,
        min_base_qual=min_base_qual,
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
    ctx: typer.Context,
    template_fasta: Path = typer.Option(
        ...,
        "--template-fasta",
        "-t",
        exists=True,
        readable=True,
        help="Template FASTA file for consensus generation.",
    ),
    config_csv: Path = typer.Option(
        ...,
        "--config-csv",
        "-C",
        exists=True,
        readable=True,
        help="CSV describing UMI/gene flanks and length bounds.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory where UMI hunter outputs will be stored.",
    ),
    umi_identity_threshold: float = typer.Option(
        0.9,
        "--umi-identity-threshold",
        "-u",
        min=0.0,
        max=1.0,
        help="UMI clustering identity threshold (default: 0.9).",
    ),
    consensus_mutation_threshold: float = typer.Option(
        0.7,
        "--consensus-mutation-threshold",
        "-M",
        min=0.0,
        max=1.0,
        help="Mutation threshold for consensus calling (default: 0.7).",
    ),
    min_cluster_size: int = typer.Option(
        1,
        "--min-cluster-size",
        "-s",
        min=1,
        help="Minimum number of reads required in a UMI cluster before a consensus is generated.",
    ),
    min_flank_ratio: int = typer.Option(
        80,
        "--min-flank-ratio",
        "-r",
        min=0,
        max=100,
        help="Minimum fuzzy match ratio (0-100) for flank detection during read extraction.",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
):
    """Cluster UMIs and generate consensus sequences from long-read FASTQ data."""
    # Validate required external tools
    try:
        validate_workflow_tools("umi_hunter")
    except ToolNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    # Apply config defaults
    config = ctx.obj.get("config", {}) if ctx.obj else {}
    umi_identity_threshold = get_option(
        config, "umi_identity_threshold", umi_identity_threshold, default=0.9, workflow="umi_hunter"
    )
    consensus_mutation_threshold = get_option(
        config, "consensus_mutation_threshold", consensus_mutation_threshold, default=0.7, workflow="umi_hunter"
    )
    min_cluster_size = get_option(
        config, "min_cluster_size", min_cluster_size, default=1, workflow="umi_hunter"
    )
    min_flank_ratio = get_option(
        config, "min_flank_ratio", min_flank_ratio, default=80, workflow="umi_hunter"
    )

    fastq_files = expand_fastq_inputs_umi(fastq)
    results = run_umi_hunter(
        template_fasta=template_fasta,
        config_csv=config_csv,
        fastq_files=fastq_files,
        output_dir=output_dir,
        umi_identity_threshold=umi_identity_threshold,
        consensus_mutation_threshold=consensus_mutation_threshold,
        min_cluster_size=min_cluster_size,
        min_flank_ratio=min_flank_ratio,
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
                f"(from {total_clusters} total) -> {entry['directory']}"
            )


@app.command("ep-library-profile", help="Profile mutation rates for ep-library sequencing data.")
def ep_library_profile_command(
    ctx: typer.Context,
    region_fasta: Path = typer.Option(
        ...,
        "--region-fasta",
        "-R",
        exists=True,
        readable=True,
        help="FASTA file describing the region of interest.",
    ),
    plasmid_fasta: Path = typer.Option(
        ...,
        "--plasmid-fasta",
        "-p",
        exists=True,
        readable=True,
        help="FASTA file with the full plasmid sequence.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help=(
            "Directory for per-sample outputs. Must be inside a workspace containing "
            "a '.uht_tooling_workspace' sentinel file."
        ),
    ),
    work_dir: Optional[Path] = typer.Option(
        None,
        "--work-dir",
        "-w",
        dir_okay=True,
        writable=True,
        help=(
            "Optional scratch directory for intermediate files (defaults to output/tmp). "
            "Must be inside a workspace containing a '.uht_tooling_workspace' sentinel file."
        ),
    ),
):
    """Quantify mutation rates for ep-library sequencing experiments."""
    # Validate required external tools
    try:
        validate_workflow_tools("ep_library_profile")
    except ToolNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

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


@app.command("ssm-profiler", help="Profile site-saturation mutagenesis libraries at target codons.")
def ssm_profiler_command(
    ctx: typer.Context,
    region_fasta: Path = typer.Option(
        ...,
        "--region-fasta",
        "-R",
        exists=True,
        readable=True,
        help="Coding-sequence FASTA for the region of interest.",
    ),
    plasmid_fasta: Path = typer.Option(
        ...,
        "--plasmid-fasta",
        "-p",
        exists=True,
        readable=True,
        help="FASTA file with the full plasmid sequence.",
    ),
    target_site: list[int] = typer.Option(
        ...,
        "--target-site",
        "-t",
        help="Target amino-acid position(s), 1-based relative to the ROI CDS. Repeat for multiple sites.",
    ),
    site_scheme: list[str] = typer.Option(
        [],
        "--site-scheme",
        "-s",
        help="Optional mapping of site to degenerate codon scheme, e.g. 45:NNK. Repeat as needed.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help=(
            "Directory for per-sample outputs. Must be inside a workspace containing "
            "a '.uht_tooling_workspace' sentinel file."
        ),
    ),
    work_dir: Optional[Path] = typer.Option(
        None,
        "--work-dir",
        "-w",
        dir_okay=True,
        writable=True,
        help=(
            "Optional scratch directory for intermediate files (defaults to output/tmp). "
            "Must be inside a workspace containing a '.uht_tooling_workspace' sentinel file."
        ),
    ),
    min_base_qual: int = typer.Option(
        0,
        "--min-base-qual",
        "-Q",
        min=0,
        max=93,
        help=(
            "Minimum per-base Phred quality score required at every position of a target "
            "codon for that read to count toward the codon read-out (0 disables filtering)."
        ),
    ),
):
    """Quantify target-site amino-acid composition in SSM libraries."""
    try:
        validate_workflow_tools("ssm_profiler")
    except ToolNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    config = ctx.obj.get("config", {}) if ctx.obj else {}
    min_base_qual = get_option(config, "min_base_qual", min_base_qual, default=0, workflow="ssm_profiler")

    fastq_files = expand_fastq_inputs_ssm(fastq)
    try:
        scheme_map = parse_scheme_map(site_scheme)
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1)

    results = run_ssm_profiler(
        fastq_paths=fastq_files,
        region_fasta=region_fasta,
        plasmid_fasta=plasmid_fasta,
        output_dir=output_dir,
        target_sites=target_site,
        scheme_map=scheme_map,
        work_dir=work_dir,
        min_base_qual=min_base_qual,
    )
    samples = results.get("samples", [])
    if not samples:
        typer.echo("No SSM profiler outputs generated.")
    else:
        typer.echo(f"Master summary written to {results['master_summary']}")
        for sample in samples:
            typer.echo(f"  Sample {sample['sample']}: {sample['results_dir']}")


@app.command("profile-inserts", help="Extract and profile inserts using probe pairs.")
def profile_inserts_command(
    ctx: typer.Context,
    probes_csv: Path = typer.Option(
        ...,
        "--probes-csv",
        "-P",
        exists=True,
        readable=True,
        help="CSV file containing upstream/downstream probes.",
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help="One or more FASTQ(.gz) paths or glob patterns (multiple --fastq options allowed).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help="Directory for per-sample outputs.",
    ),
    min_ratio: int = typer.Option(
        80,
        "--min-ratio",
        "-r",
        min=0,
        max=100,
        help="Minimum fuzzy match ratio for probe detection (default: 80).",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
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


@app.command(
    "enrich-monitor",
    help="Quantify variant enrichment from Nanopore long-read data.",
)
def enrich_monitor_command(
    ctx: typer.Context,
    reference_fasta: Path = typer.Option(
        ...,
        "--reference-fasta",
        "-R",
        exists=True,
        readable=True,
        help=(
            "Multi-sequence FASTA defining all variants to count against. "
            "Must contain at least two records."
        ),
    ),
    fastq: list[str] = typer.Option(
        ...,
        "--fastq",
        "-q",
        help=(
            "One or more FASTQ(.gz) paths or glob patterns "
            "(provide multiple --fastq options as needed)."
        ),
    ),
    baseline_sample: str = typer.Option(
        ...,
        "--baseline-sample",
        "-b",
        help=(
            "Sample name of the unselected / input-pool control. "
            "Must match the FASTQ file stem, e.g. 'control' for 'control.fastq.gz'."
        ),
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        dir_okay=True,
        writable=True,
        help=(
            "Directory for results. Must be inside a workspace containing "
            "a '.uht_tooling_workspace' sentinel file."
        ),
    ),
    target_id: Optional[str] = typer.Option(
        None,
        "--target-id",
        "-t",
        help=(
            "FASTA record ID of the target (positive-control) sequence. "
            "Defaults to the first record in the reference FASTA."
        ),
    ),
    background_id: Optional[str] = typer.Option(
        None,
        "--background-id",
        "-n",
        help=(
            "FASTA record ID of the background (negative-control) sequence. "
            "Defaults to the second record when the reference contains exactly two records."
        ),
    ),
    work_dir: Optional[Path] = typer.Option(
        None,
        "--work-dir",
        "-w",
        dir_okay=True,
        writable=True,
        help="Scratch directory for intermediate alignment files (defaults to output-dir/tmp).",
    ),
    log_path: Optional[Path] = typer.Option(
        None,
        "--log-path",
        "-l",
        dir_okay=False,
        writable=True,
        help="Optional path to write a dedicated log file.",
    ),
    n_bootstrap: int = typer.Option(
        2000,
        "--n-bootstrap",
        "-B",
        min=100,
        help=(
            "Number of bootstrap resampling iterations for Baret η and Zinchenko η' "
            "confidence intervals (default 2000)."
        ),
    ),
):
    """Quantify variant enrichment from Nanopore long-read sequencing data."""
    try:
        validate_workflow_tools("enrich_monitor")
    except ToolNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    fastq_files = expand_fastq_inputs_enrich(fastq)
    try:
        results = run_enrich_monitor(
            reference_fasta=reference_fasta,
            fastq_paths=fastq_files,
            baseline_sample=baseline_sample,
            output_dir=output_dir,
            target_id=target_id,
            background_id=background_id,
            work_dir=work_dir,
            log_path=log_path,
            n_bootstrap=n_bootstrap,
        )
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1)

    typer.echo(f"Enrichment stats: {results['enrichment_stats']}")
    typer.echo(f"Abundance table:  {results['abundance_table']}")
    typer.echo(f"Summary figure:   {results['summary_figure']}")
    for entry in results["samples"]:
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
    legacy: bool = typer.Option(
        False,
        "--legacy",
        help="Use the legacy Gradio GUI instead of the new NiceGUI frontend.",
    ),
    share: bool = typer.Option(
        False,
        "--share",
        help="Enable Gradio's public sharing tunnel (legacy mode only).",
    ),
):
    """Launch the graphical interface (NiceGUI by default, --legacy for Gradio)."""
    try:
        if legacy:
            from uht_tooling.workflows.gui import launch_gui
            launch_gui(server_name=server_name, server_port=server_port, share=share)
        else:
            launch_web_gui(host=server_name, port=server_port)
    except KeyboardInterrupt:
        typer.echo("GUI stopped by user.")
    except Exception as exc:
        typer.echo(f"Failed to start GUI: {exc}")
        raise typer.Exit(1)


def main():
    app()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Gradio GUI for the uht_tooling package built on the refactored workflows."""

import contextlib
import logging
import os
import socket
import shutil
import tempfile
import textwrap
import zipfile
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple

try:
    import gradio as gr
    import pandas as pd
except ImportError as exc:  # pragma: no cover - handled at runtime
    raise SystemExit(
        "Missing dependency: "
        f"{exc}. Install optional GUI extras via 'pip install gradio pandas'."
    ) from exc

from uht_tooling.workflows.design_gibson import run_design_gibson
from uht_tooling.workflows.design_slim import run_design_slim
from uht_tooling.workflows.mut_rate import run_ep_library_profile
from uht_tooling.workflows.mutation_caller import run_mutation_caller
from uht_tooling.workflows.nextera_designer import run_nextera_primer_design
from uht_tooling.workflows.profile_inserts import run_profile_inserts
from uht_tooling.workflows.umi_hunter import run_umi_hunter

_LOGGER = logging.getLogger("uht_tooling.gui")


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _ensure_text(value: str, field: str) -> str:
    value = (value or "").strip()
    if not value:
        raise ValueError(f"{field} cannot be empty.")
    return value


def _clean_temp_path(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path, ignore_errors=True)


def _zip_paths(paths: Iterable[Path], prefix: str) -> Path:
    archive_dir = Path(tempfile.mkdtemp(prefix=f"uht_gui_{prefix}_zip_"))
    zip_path = archive_dir / f"{prefix}_results.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            path = Path(path)
            if not path.exists():
                continue
            if path.is_dir():
                for file in path.rglob("*"):
                    if file.is_file():
                        arcname = Path(path.name) / file.relative_to(path)
                        archive.write(file, arcname.as_posix())
            else:
                archive.write(path, arcname=path.name)
    return zip_path


def _preview_csv(csv_path: Path, max_rows: int = 10) -> str:
    if not csv_path.exists():
        return "*(output file missing)*"
    try:
        df = pd.read_csv(csv_path)
    except Exception as exc:  # pragma: no cover - defensive
        return f"*(unable to read CSV: {exc})*"
    if df.empty:
        return "*(no rows generated)*"
    return df.head(max_rows).to_markdown(index=False)


def _format_header(title: str) -> str:
    return f"### {title}\n"


def _port_is_available(host: str, port: int) -> bool:
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind((host, port))
        except OSError:
            return False
    return True


def _find_server_port(host: str, preferred: Optional[int]) -> int:
    env_port = os.getenv("GRADIO_SERVER_PORT")
    if env_port:
        try:
            env_value = int(env_port)
        except ValueError:
            _LOGGER.warning("Invalid GRADIO_SERVER_PORT=%s; ignoring.", env_port)
        else:
            preferred = env_value

    if preferred is not None and _port_is_available(host, preferred):
        return preferred

    if preferred is not None:
        _LOGGER.warning(
            "Preferred port %s is unavailable on %s. Searching for an open port.",
            preferred,
            host,
        )

    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind((host, 0))
        return sock.getsockname()[1]


# ---------------------------------------------------------------------------
# Workflow adapters used by the GUI tabs
# ---------------------------------------------------------------------------

def run_gui_nextera(forward_primer: str, reverse_primer: str) -> Tuple[str, Optional[str]]:
    try:
        forward = _ensure_text(forward_primer, "Forward primer")
        reverse = _ensure_text(reverse_primer, "Reverse primer")

        work_dir = Path(tempfile.mkdtemp(prefix="uht_gui_nextera_work_"))
        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_nextera_out_"))

        binding_csv = work_dir / "binding_regions.csv"
        binding_csv.write_text("binding_region\n" + forward + "\n" + reverse + "\n")

        output_csv = output_dir / "nextera_xt_primers.csv"
        result_csv = run_nextera_primer_design(binding_csv, output_csv)

        summary = _format_header("Nextera XT Primers") + _preview_csv(result_csv)
        archive = _zip_paths([output_dir], "nextera")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover - runtime feedback
        _LOGGER.exception("Nextera GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        _clean_temp_path(locals().get("work_dir", Path()))
        _clean_temp_path(locals().get("output_dir", Path()))


def run_gui_design_slim(
    template_gene_content: str,
    context_content: str,
    mutations_text: str,
) -> Tuple[str, Optional[str]]:
    try:
        gene_seq = _ensure_text(template_gene_content, "Template gene sequence")
        context_seq = _ensure_text(context_content, "Context sequence")
        mutation_lines = [line.strip() for line in mutations_text.splitlines() if line.strip()]
        if not mutation_lines:
            raise ValueError("Provide at least one mutation (e.g., A123G).")

        work_dir = Path(tempfile.mkdtemp(prefix="uht_gui_slim_work_"))
        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_slim_out_"))

        gene_fasta = work_dir / "template_gene.fasta"
        context_fasta = work_dir / "context.fasta"
        mutations_csv = work_dir / "mutations.csv"

        gene_fasta.write_text(f">template\n{gene_seq}\n")
        context_fasta.write_text(f">context\n{context_seq}\n")
        mutations_csv.write_text("mutations\n" + "\n".join(mutation_lines) + "\n")

        result_csv = run_design_slim(
            gene_fasta=gene_fasta,
            context_fasta=context_fasta,
            mutations_csv=mutations_csv,
            output_dir=output_dir,
        )

        summary = _format_header("SLIM Primer Design") + _preview_csv(result_csv)
        archive = _zip_paths([output_dir], "slim")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("SLIM GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        _clean_temp_path(locals().get("work_dir", Path()))
        _clean_temp_path(locals().get("output_dir", Path()))


def run_gui_design_gibson(
    template_gene_content: str,
    context_content: str,
    mutations_text: str,
) -> Tuple[str, Optional[str]]:
    try:
        gene_seq = _ensure_text(template_gene_content, "Template gene sequence")
        context_seq = _ensure_text(context_content, "Context sequence")
        mutation_lines = [line.strip() for line in mutations_text.splitlines() if line.strip()]
        if not mutation_lines:
            raise ValueError("Provide at least one mutation (e.g., A123G).")

        work_dir = Path(tempfile.mkdtemp(prefix="uht_gui_gibson_work_"))
        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_gibson_out_"))

        gene_fasta = work_dir / "template_gene.fasta"
        context_fasta = work_dir / "context.fasta"
        mutations_csv = work_dir / "mutations.csv"

        gene_fasta.write_text(f">template\n{gene_seq}\n")
        context_fasta.write_text(f">context\n{context_seq}\n")
        mutations_csv.write_text("mutations\n" + "\n".join(mutation_lines) + "\n")

        outputs = run_design_gibson(
            gene_fasta=gene_fasta,
            context_fasta=context_fasta,
            mutations_csv=mutations_csv,
            output_dir=output_dir,
        )

        primers_csv = Path(outputs["primers_csv"])
        plan_csv = Path(outputs["plan_csv"])
        summary = _format_header("Gibson Assembly") + "\n".join(
            [
                "**Primer preview**",
                _preview_csv(primers_csv),
                "",
                "**Assembly plan preview**",
                _preview_csv(plan_csv),
            ]
        )
        archive = _zip_paths([output_dir], "gibson")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("Gibson GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        _clean_temp_path(locals().get("work_dir", Path()))
        _clean_temp_path(locals().get("output_dir", Path()))


def run_gui_mutation_caller(
    fastq_file: Optional[str],
    template_file: Optional[str],
    upstream_flank: str,
    downstream_flank: str,
    min_gene_length: Optional[float],
    max_gene_length: Optional[float],
) -> Tuple[str, Optional[str]]:
    config_dir: Optional[Path] = None
    output_dir: Optional[Path] = None
    try:
        if not fastq_file or not template_file:
            raise ValueError("Upload a FASTQ(.gz) read file and the reference template FASTA.")

        gene_start = _ensure_text(upstream_flank, "Upstream flank")
        gene_end = _ensure_text(downstream_flank, "Downstream flank")
        if min_gene_length is None or max_gene_length is None:
            raise ValueError("Provide minimum and maximum gene lengths (in nucleotides).")

        gene_min = int(min_gene_length)
        gene_max = int(max_gene_length)
        if gene_min <= 0 or gene_max <= 0:
            raise ValueError("Gene length bounds must be positive integers.")
        if gene_min > gene_max:
            raise ValueError("Minimum gene length cannot exceed the maximum gene length.")

        config_dir = Path(tempfile.mkdtemp(prefix="uht_gui_mutation_cfg_"))
        config_csv = config_dir / "mutation_flanks.csv"
        pd.DataFrame(
            {
                "gene_flanks": [gene_start.upper(), gene_end.upper()],
                "gene_min_max": [gene_min, gene_max],
            }
        ).to_csv(config_csv, index=False)

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_mutation_out_"))
        results = run_mutation_caller(
            template_fasta=Path(template_file),
            flanks_csv=config_csv,
            fastq_files=[Path(fastq_file)],
            output_dir=output_dir,
            threshold=10,
        )

        if not results:
            return "No amino-acid substitutions detected. Check flank selections and read quality.", None

        lines = [
            "### Mutation Caller",
            "",
            "Long-read reads were aligned to the provided template, flank-delimited coding regions were extracted, and amino-acid substitutions were summarised.",
            "",
            "**Run outputs**",
        ]
        sample_dirs = []
        for entry in results:
            lines.append(f"- **{entry['sample']}** → {entry['directory']}")
            sample_dirs.append(Path(entry["directory"]))
        summary = "\n".join(lines)
        archive = _zip_paths(sample_dirs, "mutation_caller")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("Mutation caller GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        if config_dir:
            _clean_temp_path(config_dir)
        if output_dir:
            _clean_temp_path(output_dir)


def run_gui_umi_hunter(
    fastq_file: Optional[str],
    template_file: Optional[str],
    umi_start: str,
    umi_end: str,
    umi_min_length: Optional[float],
    umi_max_length: Optional[float],
    gene_start: str,
    gene_end: str,
    umi_identity_threshold: float,
    consensus_threshold: float,
    min_cluster_size: int,
) -> Tuple[str, Optional[str]]:
    config_dir: Optional[Path] = None
    output_dir: Optional[Path] = None
    try:
        if not fastq_file or not template_file:
            raise ValueError("Upload a FASTQ(.gz) read file and the template FASTA.")

        umi_start_clean = _ensure_text(umi_start, "UMI upstream flank").upper()
        umi_end_clean = _ensure_text(umi_end, "UMI downstream flank").upper()
        gene_start_clean = _ensure_text(gene_start, "Gene upstream flank").upper()
        gene_end_clean = _ensure_text(gene_end, "Gene downstream flank").upper()
        if umi_min_length is None or umi_max_length is None:
            raise ValueError("Provide minimum and maximum UMI lengths.")

        umi_min = int(umi_min_length)
        umi_max = int(umi_max_length)
        if umi_min <= 0 or umi_max <= 0:
            raise ValueError("UMI length bounds must be positive integers.")
        if umi_min > umi_max:
            raise ValueError("Minimum UMI length cannot exceed the maximum length.")
        if not (0.0 <= umi_identity_threshold <= 1.0):
            raise ValueError("UMI identity threshold must be between 0 and 1.")
        if not (0.0 <= consensus_threshold <= 1.0):
            raise ValueError("Consensus mutation threshold must be between 0 and 1.")
        if min_cluster_size is None or int(min_cluster_size) < 1:
            raise ValueError("Minimum cluster size must be at least 1.")
        min_cluster_size_int = int(min_cluster_size)

        config_dir = Path(tempfile.mkdtemp(prefix="uht_gui_umi_cfg_"))
        config_csv = config_dir / "umi_config.csv"
        pd.DataFrame(
            {
                "umi_flanks": [umi_start_clean, umi_end_clean],
                "umi_min_max": [umi_min, umi_max],
                "gene_flanks": [gene_start_clean, gene_end_clean],
            }
        ).to_csv(config_csv, index=False)

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_umi_out_"))
        results = run_umi_hunter(
            template_fasta=Path(template_file),
            config_csv=config_csv,
            fastq_files=[Path(fastq_file)],
            output_dir=output_dir,
            umi_identity_threshold=umi_identity_threshold,
            consensus_mutation_threshold=consensus_threshold,
            min_cluster_size=min_cluster_size_int,
        )

        if not results:
            return (
                "No UMI clusters were generated. Double-check flank selections and threshold settings.",
                None,
            )

        lines = [
            "### UMI Hunter",
            "",
            "Reads were scanned for UMI and gene flanks, deduplicated by UMI, and consensus alleles were generated.",
            "",
            "**Run outputs**",
        ]
        sample_dirs = []
        for entry in results:
            total_clusters = entry.get("clusters_total", entry["clusters"])
            lines.append(
                f"- **{entry['sample']}** → {entry['clusters']} consensus clusters "
                f"(≥ {min_cluster_size_int} reads) from {total_clusters} total, "
                f"results in {entry['directory']}"
            )
            sample_dirs.append(Path(entry["directory"]))
        summary = "\n".join(lines)
        archive = _zip_paths(sample_dirs, "umi_hunter")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("UMI hunter GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        if config_dir:
            _clean_temp_path(config_dir)
        if output_dir:
            _clean_temp_path(output_dir)


def run_gui_profile_inserts(
    probes_table: Any,
    fastq_files: Sequence[str],
    min_ratio: int,
) -> Tuple[str, Optional[str]]:
    config_dir: Optional[Path] = None
    output_dir: Optional[Path] = None
    try:
        if not fastq_files:
            raise ValueError("Upload at least one FASTQ(.gz) file.")
        if probes_table is None:
            raise ValueError("Provide at least one probe pair.")

        if isinstance(probes_table, pd.DataFrame):
            df = probes_table.copy()
        else:
            df = pd.DataFrame(probes_table or [], columns=["name", "upstream", "downstream"])

        # Normalise and validate probe entries
        df = df.replace({pd.NA: "", None: ""})
        for column in df.columns:
            if df[column].dtype == object:
                df[column] = df[column].map(lambda x: x.strip() if isinstance(x, str) else x)

        if "upstream" not in df.columns or "downstream" not in df.columns:
            raise ValueError("Probe table must contain 'upstream' and 'downstream' columns.")

        df_valid = df[(df["upstream"] != "") & (df["downstream"] != "")].copy()
        if df_valid.empty:
            raise ValueError("Enter at least one probe pair with both upstream and downstream sequences.")

        df_valid = df_valid.reset_index(drop=True)
        if "name" not in df_valid.columns:
            df_valid["name"] = [f"probe_{i + 1}" for i in range(len(df_valid))]
        else:
            fallback_names = pd.Series(
                [f"probe_{i + 1}" for i in range(len(df_valid))], index=df_valid.index
            )
            df_valid["name"] = df_valid["name"].replace("", pd.NA).fillna(fallback_names)

        config_dir = Path(tempfile.mkdtemp(prefix="uht_gui_profile_cfg_"))
        probes_csv = config_dir / "probes.csv"
        df_valid.to_csv(probes_csv, index=False)

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_profile_out_"))
        results = run_profile_inserts(
            probes_csv=probes_csv,
            fastq_files=[Path(f) for f in fastq_files],
            output_dir=output_dir,
            min_ratio=int(min_ratio),
        )

        if not results:
            return "No inserts were extracted. Adjust probe sequences or similarity threshold and try again.", None

        first_insert = results[0]["fasta"] if isinstance(results, list) else None
        preview = "*(preview unavailable)*"
        if first_insert and Path(first_insert).exists():
            preview = Path(first_insert).read_text().splitlines()[0][:120] + "..."

        summary = textwrap.dedent(
            """
            ### Insert Profiling
            Probe-defined regions were scanned in the provided FASTQ files, inserts were extracted, and QC metrics were generated.

            **Key outputs**
            - FASTA files containing extracted inserts per probe pair
            - Summary tables covering length, GC content, duplicate rate, and probe match quality
            - A gallery of QC plots (length distributions, base composition, probe performance)
            """
        )
        archive = _zip_paths([Path(r["directory"]) for r in results], "profile_inserts")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("Profile inserts GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        if config_dir:
            _clean_temp_path(config_dir)
        if output_dir:
            _clean_temp_path(output_dir)


def run_gui_ep_library_profile(
    fastq_files: Sequence[str],
    region_fasta: Optional[str],
    plasmid_fasta: Optional[str],
) -> Tuple[str, Optional[str]]:
    try:
        if not fastq_files or not region_fasta or not plasmid_fasta:
            raise ValueError("Upload FASTQ(.gz) files plus region-of-interest and plasmid FASTA files.")

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_ep_out_"))
        work_dir = Path(tempfile.mkdtemp(prefix="uht_gui_ep_work_"))
        results = run_ep_library_profile(
            fastq_paths=[Path(f) for f in fastq_files],
            region_fasta=Path(region_fasta),
            plasmid_fasta=Path(plasmid_fasta),
            output_dir=output_dir,
            work_dir=work_dir,
        )

        master_summary = Path(results["master_summary"])
        summary_text = master_summary.read_text() if master_summary.exists() else "Summary unavailable."

        lines = ["### EP Library Profile", "", "**Master summary**", "```", summary_text.strip(), "```"]
        for sample in results.get("samples", []):
            lines.append(f"- {sample['sample']} → {sample['results_dir']}")
        summary = "\n".join(lines)

        sample_dirs = [Path(sample["results_dir"]) for sample in results.get("samples", [])]
        archive = _zip_paths(sample_dirs + [master_summary], "ep_library")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("EP library profile GUI failure")
        return f"⚠️ Error: {exc}", None
    finally:
        _clean_temp_path(locals().get("work_dir", Path()))


# ---------------------------------------------------------------------------
# Gradio Interface
# ---------------------------------------------------------------------------

def create_gui() -> gr.Blocks:
    custom_css = """
    .gradio-container {
        max-width: 1400px;
        margin: 0 auto;
    }
    .hero {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 12px;
        color: white;
        margin-bottom: 2rem;
        text-align: center;
    }
    """

    with gr.Blocks(title="uht-tooling GUI", css=custom_css) as demo:
        with gr.Column(elem_classes="hero"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    # uht-tooling
                    A guided graphical interface for primer design and sequencing analysis. Each tab mirrors the command-line workflows documented in the README and bundles results, logs, and QC artefacts for download.

                    **How to use**
                    1. Select the workflow that matches your experiment.
                    2. Provide the required inputs (text fields, FASTQ/FASTA uploads, or probe tables).
                    3. Run the analysis and download the ZIP archive for complete outputs.

                    Need automation or batch processing? Use the Typer CLI (`uht-tooling ...`) with the same arguments shown here.
                    """
                )
            )

        with gr.Tab("Nextera XT"):  # --- Nextera ---
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Illumina-Compatible Primer Design
                    Generates Nextera XT-ready primers from forward/reverse binding regions. The workflow preloads 12 i5 and 12 i7 indices (144 combinations) and mirrors the “One-PCR-to-flowcell” process described in the README.

                    **Inputs**
                    - Forward primer binding region (5'→3')
                    - Reverse primer binding region (5'→3')

                    **Outputs**
                    - CSV with i5/i7 indices, primer sequences, and ordering-ready metadata.
                    - Run log noting index selection and any validation warnings.
                    """
                )
            )
            forward = gr.Textbox(label="Forward primer (5'→3')")
            reverse = gr.Textbox(label="Reverse primer (5'→3')")
            nextera_btn = gr.Button("Generate Primers", variant="primary")
            nextera_summary = gr.Markdown(label="Summary")
            nextera_download = gr.File(label="Download primers", file_count="single")
            nextera_btn.click(
                fn=run_gui_nextera,
                inputs=[forward, reverse],
                outputs=[nextera_summary, nextera_download],
            )
            with gr.Accordion("Wet-lab guidance", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - Monitor amplification by qPCR and cap the cycle count to reach roughly 10 % yield to limit bias.
                        - Purify products with SPRIselect beads (~0.65:1 bead:DNA ratio) to remove residual primers.
                        - Confirm primer depletion via electrophoresis (e.g., BioAnalyzer) before sequencing prep.
                        """
                    )
                )

        with gr.Tab("SLIM"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Sequence-Ligation Independent Mutagenesis
                    Designs paired short/long primers to introduce targeted mutations by SLIM cloning, matching the workflow outlined in the README.

                    **Inputs**
                    - Target gene coding sequence (FASTA content).
                    - Plasmid or genomic context containing the gene.
                    - Mutations (one per line, e.g. substitution `A123G`, deletion `T241Del`, insertion `T241TS`).

                    **Outputs**
                    - `SLIM_primers.csv` with primer sequences and annealing temperatures.
                    - Log file capturing primer QC and any design warnings.
                    """
                )
            )
            slim_gene = gr.Textbox(label="Gene sequence", lines=4)
            slim_context = gr.Textbox(label="Plasmid context", lines=4)
            slim_mutations = gr.Textbox(label="Mutations (one per line)", lines=6)
            slim_btn = gr.Button("Design SLIM primers", variant="primary")
            slim_summary = gr.Markdown(label="Summary")
            slim_download = gr.File(label="Download primers", file_count="single")
            slim_btn.click(
                fn=run_gui_design_slim,
                inputs=[slim_gene, slim_context, slim_mutations],
                outputs=[slim_summary, slim_download],
            )
            with gr.Accordion("Bench workflow blueprint", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        1. Run two PCRs: (A) long forward + short reverse, (B) long reverse + short forward.
                        2. Combine 10 µL from each PCR with 10 µL H-buffer (150 mM Tris pH 8, 400 mM NaCl, 60 mM EDTA).
                        3. Thermocycle: 99 °C 3 min → 2× (65 °C 5 min → 30 °C 15 min) → hold at 4 °C.
                        4. Transform directly into NEB 5-alpha or BL21 (DE3); the method scales to dozens of mutants simultaneously.
                        """
                    )
                )

        with gr.Tab("Gibson"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Gibson Assembly Primer Design
                    Plans primer sets and assembly steps for Gibson mutagenesis, supporting multi-mutation constructs using the `+` syntax (e.g. `A123G+T150A`).

                    **Inputs**
                    - Coding sequence for the gene of interest.
                    - Circular plasmid context sequence.
                    - Mutation definitions (one per line; use `+` to bundle simultaneous edits).

                    **Outputs**
                    - Primer CSV with overlap sequences and melting temperatures.
                    - Assembly plan CSV detailing fragment combinations.
                    - Log summarising design decisions and any warnings about overlapping regions.
                    """
                )
            )
            gibson_gene = gr.Textbox(label="Gene sequence", lines=4)
            gibson_context = gr.Textbox(label="Plasmid context", lines=4)
            gibson_mutations = gr.Textbox(label="Mutations", lines=6)
            gibson_btn = gr.Button("Design Gibson primers", variant="primary")
            gibson_summary = gr.Markdown(label="Summary")
            gibson_download = gr.File(label="Download results", file_count="single")
            gibson_btn.click(
                fn=run_gui_design_gibson,
                inputs=[gibson_gene, gibson_context, gibson_mutations],
                outputs=[gibson_summary, gibson_download],
            )
            with gr.Accordion("Tips for multi-mutation designs", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - If two mutations compete for primer space, design them in sequential runs to avoid overly long primers.
                        - Use the assembly plan CSV to map which fragments to combine in each Gibson reaction.
                        - When replacing entire codons (e.g. `L46GP`), ensure the plasmid context covers both flanks to maintain overlap.
                        """
                    )
                )

        with gr.Tab("Mutation Caller"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Long-read Mutation Analysis
                    Extracts coding regions bounded by user-defined flanks, aligns them to the template, and reports amino-acid substitutions alongside co-occurrence summaries.

                    **Required inputs**
                    - FASTQ (.fastq.gz): Oxford Nanopore or other long-read data.
                    - Template FASTA: coding sequence used as the reference for alignment.
                    - Flank sequences: short 8–12 bp motifs immediately upstream and downstream of the gene.
                    - Gene length bounds: acceptable size window (in nucleotides) for the extracted gene segment.
                    """
                )
            )
            with gr.Row():
                mc_fastq = gr.File(
                    label="FASTQ (.fastq.gz)",
                    file_types=[".fastq", ".gz"],
                    type="filepath",
                )
                mc_template = gr.File(
                    label="Template FASTA",
                    file_types=[".fasta", ".fa"],
                    type="filepath",
                )
            with gr.Row():
                mc_upstream = gr.Textbox(
                    label="Upstream flank (5'→3')",
                    placeholder="e.g. ACTGTTAG",
                )
                mc_downstream = gr.Textbox(
                    label="Downstream flank (5'→3')",
                    placeholder="e.g. CGAACCTA",
                )
            with gr.Row():
                mc_min_len = gr.Number(
                    label="Minimum gene length (nt)",
                    value=900,
                    precision=0,
                )
                mc_max_len = gr.Number(
                    label="Maximum gene length (nt)",
                    value=1200,
                    precision=0,
                )
            mc_btn = gr.Button("Run mutation caller", variant="primary")
            mc_summary = gr.Markdown(label="Summary")
            mc_download = gr.File(label="Download results", file_count="single")
            mc_btn.click(
                fn=run_gui_mutation_caller,
                inputs=[
                    mc_fastq,
                    mc_template,
                    mc_upstream,
                    mc_downstream,
                    mc_min_len,
                    mc_max_len,
                ],
                outputs=[mc_summary, mc_download],
            )
            with gr.Accordion("What happens under the hood", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - Reads are scanned for the upstream and downstream flanks; the sequence between them is treated as the gene of interest if it falls within the specified length window.
                        - MAFFT aligns recovered genes to the reference template and the pipeline annotates amino-acid substitutions, co-occurrence networks, and depth statistics.
                        - Outputs mirror the CLI version: per-sample directories with CSV summaries, JSON co-occurrence graphs, QC plots, and a detailed `run.log`.
                        """
                    )
                )

        with gr.Tab("UMI Hunter"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### UMI–Gene Pair Clustering
                    Detects UMI barcodes, extracts paired gene inserts, clusters reads by UMI identity, and emits consensus sequences with abundance tables.

                    **Required inputs**
                    - FASTQ (.fastq.gz) containing UMI-tagged reads.
                    - Template FASTA for downstream consensus calling.
                    - UMI and gene flank sequences marking the barcode and insert boundaries.
                    - UMI length bounds plus clustering thresholds.
                    - Minimum reads per cluster to keep (clusters below the threshold are reported but no consensus is generated).
                    """
                )
            )
            with gr.Row():
                umi_fastq = gr.File(
                    label="FASTQ (.fastq.gz)",
                    file_types=[".fastq", ".gz"],
                    type="filepath",
                )
                umi_template = gr.File(
                    label="Template FASTA",
                    file_types=[".fasta", ".fa"],
                    type="filepath",
                )
            with gr.Row():
                umi_start = gr.Textbox(
                    label="UMI upstream flank (5'→3')",
                    placeholder="e.g. ACACTCTTTCCCTACACGAC",
                )
                umi_end = gr.Textbox(
                    label="UMI downstream flank (5'→3')",
                    placeholder="e.g. GACTGGAGTTCAGACGTGTG",
                )
            with gr.Row():
                gene_start = gr.Textbox(
                    label="Gene upstream flank (5'→3')",
                    placeholder="e.g. ATG...",
                )
                gene_end = gr.Textbox(
                    label="Gene downstream flank (5'→3')",
                    placeholder="e.g. TTA...",
                )
            with gr.Row():
                umi_min_len = gr.Number(
                    label="Minimum UMI length (nt)",
                    value=8,
                    precision=0,
                )
                umi_max_len = gr.Number(
                    label="Maximum UMI length (nt)",
                    value=14,
                    precision=0,
                )
            with gr.Row():
                umi_identity = gr.Slider(
                    label="UMI clustering identity",
                    minimum=0.5,
                    maximum=1.0,
                    value=0.9,
                    step=0.05,
                )
                consensus_threshold = gr.Slider(
                    label="Consensus mutation threshold",
                    minimum=0.5,
                    maximum=1.0,
                    value=0.7,
                    step=0.05,
                )
            umi_min_cluster = gr.Slider(
                label="Minimum reads per cluster",
                minimum=1,
                maximum=50,
                value=3,
                step=1,
            )
            umi_btn = gr.Button("Run UMI hunter", variant="primary")
            umi_summary = gr.Markdown(label="Summary")
            umi_download = gr.File(label="Download results", file_count="single")
            umi_btn.click(
                fn=run_gui_umi_hunter,
                inputs=[
                    umi_fastq,
                    umi_template,
                    umi_start,
                    umi_end,
                    umi_min_len,
                    umi_max_len,
                    gene_start,
                    gene_end,
                    umi_identity,
                    consensus_threshold,
                    umi_min_cluster,
                ],
                outputs=[umi_summary, umi_download],
            )
            with gr.Accordion("What the pipeline generates", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - Reads are searched for the UMI barcode and gene flanks on both strands; valid pairs feed into UMI grouping.
                        - UMIs within the chosen identity threshold are merged, and consensus sequences are computed with the mutation threshold.
                        - Outputs include per-sample summaries, consensus FASTA files, cluster membership tables, QC plots, and logs mirroring the CLI workflow.
                        """
                    )
                )

        with gr.Tab("Profile Inserts"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Probe-Guided Insert Profiling
                    Characterises inserts demarcated by user-supplied upstream/downstream probes, extracts sequences, and produces QC plots plus summary tables.

                    **Required inputs**
                    - FASTQ reads containing the inserts of interest.
                    - One or more probe pairs: 5'→3' sequences for the upstream and downstream anchors (reverse complements are matched automatically).
                    """
                )
            )
            probes_table = gr.Dataframe(
                headers=["name (optional)", "upstream", "downstream"],
                datatype=["str", "str", "str"],
                row_count=(1, "dynamic"),
                col_count=3,
                value=[["probe_1", "", ""]],
                interactive=True,
                label="Probe pairs",
            )
            pi_fastq = gr.File(
                label="FASTQ files (.fastq/.gz)",
                file_types=[".fastq", ".gz"],
                file_count="multiple",
                type="filepath",
            )
            pi_ratio = gr.Slider(
                label="Minimum fuzzy-match ratio",
                minimum=50,
                maximum=100,
                value=80,
                step=1,
            )
            pi_btn = gr.Button("Profile inserts", variant="primary")
            pi_summary = gr.Markdown(label="Summary")
            pi_download = gr.File(label="Download results", file_count="single")
            pi_btn.click(
                fn=run_gui_profile_inserts,
                inputs=[probes_table, pi_fastq, pi_ratio],
                outputs=[pi_summary, pi_download],
            )
            with gr.Accordion("Output overview", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - Inserts are extracted whenever probe matches are detected above the chosen similarity threshold (default 80).
                        - A FASTA file of inserts, probe-level QC metrics, base composition summaries, and a suite of plots (length distribution, GC content, duplicate rate, probe performance) are packaged for each input FASTQ.
                        - Logs are stored alongside the results so runs remain fully reproducible.
                        """
                    )
                )

        with gr.Tab("EP Library Profile"):
            gr.Markdown(
                textwrap.dedent(
                    """
                    ### Library Profiling Without UMIs
                    Estimates background and target mutation rates for enzyme evolution libraries without UMI barcodes.

                    **Inputs**
                    - FASTQ reads (*.fastq/.gz) from the ep-library experiment.
                    - Region-of-interest FASTA delineating the mutational window.
                    - Plasmid FASTA providing the full reference context.

                    **Outputs**
                    - Per-sample directories with coverage tables, mutation rate statistics, and QC plots.
                    - `master_summary.txt` aggregating condition-level metrics.
                    - Verbose logs recording alignment commands and rate calculations.
                    """
                )
            )
            ep_fastq = gr.File(
                label="FASTQ files",
                file_types=[".fastq", ".gz"],
                file_count="multiple",
                type="filepath",
            )
            ep_region = gr.File(label="Region FASTA", file_types=[".fasta", ".fa"], type="filepath")
            ep_plasmid = gr.File(label="Plasmid FASTA", file_types=[".fasta", ".fa"], type="filepath")
            ep_btn = gr.Button("Run profiling", variant="primary")
            ep_summary = gr.Markdown(label="Summary")
            ep_download = gr.File(label="Download results", file_count="single")
            ep_btn.click(
                fn=run_gui_ep_library_profile,
                inputs=[ep_fastq, ep_region, ep_plasmid],
                outputs=[ep_summary, ep_download],
            )
            with gr.Accordion("How mutation rates are derived", open=False):
                gr.Markdown(
                    textwrap.dedent(
                        """
                        - Reads are aligned against both the region-of-interest and the full plasmid to measure target and background mismatch rates; their difference yields the net nucleotide mutation rate with propagated binomial and quality-score uncertainty.
                        - The net per-base rate is multiplied by the CDS length to obtain λ₍bp₎ (mutations per copy), then Monte Carlo simulations flip random bases, translate the mutated CDS, and count amino-acid differences—those simulated means and confidence intervals are the values plotted in the QC figure.
                        - When multiple Q-score thresholds are analysed, the CLI combines them via a precision-weighted consensus (after discarding filters with <1000 mappable bases). The consensus AA mutation rate is written to `aa_mutation_consensus.txt` and drawn as a horizontal guide in the plot.
                        - Download the archive to inspect per-sample plots, TSV summaries, the consensus summary, and logs for troubleshooting.
                        """
                    )
                )

        gr.Markdown(
            textwrap.dedent(
                """
                ---
                **Tips for new users**

                1. Prepare your inputs (FASTA/CSV/FASTQ) before opening the tab.
                2. Click the action button and wait for the summary to appear.
                3. Download the ZIP archive for the complete result set.
                4. For automation or batch processing, use the command-line interface instead (`uht-tooling ...`).
                """
            )
        )

    return demo


def launch_gui(
    server_name: str = "127.0.0.1",
    server_port: Optional[int] = 7860,
    share: bool = False,
) -> None:
    resolved_port = _find_server_port(server_name, server_port)
    _LOGGER.info("Starting uht-tooling GUI on http://%s:%s", server_name, resolved_port)
    demo = create_gui()
    demo.launch(
        server_name=server_name,
        server_port=resolved_port,
        share=share,
        show_error=True,
    )


def main() -> None:  # pragma: no cover - entry point wrapper
    logging.basicConfig(level=logging.INFO)
    launch_gui()


if __name__ == "__main__":  # pragma: no cover
    main()

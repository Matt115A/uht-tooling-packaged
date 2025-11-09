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
from typing import Iterable, List, Optional, Sequence, Tuple

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
    config_csv_file: Optional[str],
) -> Tuple[str, Optional[str]]:
    try:
        if not fastq_file or not template_file or not config_csv_file:
            raise ValueError("Upload a FASTQ(.gz), template FASTA, and configuration CSV.")

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_mutation_out_"))
        results = run_mutation_caller(
            template_fasta=Path(template_file),
            flanks_csv=Path(config_csv_file),
            fastq_files=[Path(fastq_file)],
            output_dir=output_dir,
            threshold=10,
        )

        if not results:
            return "No amino-acid substitutions detected.", None

        lines = ["### Mutation Caller", ""]
        sample_dirs = []
        for entry in results:
            lines.append(f"**{entry['sample']}** → {entry['directory']}")
            sample_dirs.append(Path(entry["directory"]))
        summary = "\n".join(lines)
        archive = _zip_paths(sample_dirs, "mutation_caller")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("Mutation caller GUI failure")
        return f"⚠️ Error: {exc}", None


def run_gui_umi_hunter(
    fastq_file: Optional[str],
    template_file: Optional[str],
    config_csv_file: Optional[str],
) -> Tuple[str, Optional[str]]:
    try:
        if not fastq_file or not template_file or not config_csv_file:
            raise ValueError("Upload a FASTQ(.gz), template FASTA, and configuration CSV.")

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_umi_out_"))
        results = run_umi_hunter(
            template_fasta=Path(template_file),
            config_csv=Path(config_csv_file),
            fastq_files=[Path(fastq_file)],
            output_dir=output_dir,
        )

        if not results:
            return "No UMI clusters were generated. Check input quality and thresholds.", None

        lines = ["### UMI Hunter", ""]
        sample_dirs = []
        for entry in results:
            lines.append(
                f"**{entry['sample']}** → {entry['clusters']} clusters, results in {entry['directory']}"
            )
            sample_dirs.append(Path(entry["directory"]))
        summary = "\n".join(lines)
        archive = _zip_paths(sample_dirs, "umi_hunter")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("UMI hunter GUI failure")
        return f"⚠️ Error: {exc}", None


def run_gui_profile_inserts(
    probes_csv_path: Optional[str],
    fastq_files: Sequence[str],
) -> Tuple[str, Optional[str]]:
    try:
        if not probes_csv_path or not fastq_files:
            raise ValueError("Upload the probe CSV and at least one FASTQ(.gz) file.")

        output_dir = Path(tempfile.mkdtemp(prefix="uht_gui_profile_out_"))
        results = run_profile_inserts(
            probes_csv=Path(probes_csv_path),
            fastq_files=[Path(f) for f in fastq_files],
            output_dir=output_dir,
        )

        if not results:
            return "No inserts were extracted. Adjust probe settings and try again.", None

        first_insert = results[0]["fasta"] if isinstance(results, list) else None
        preview = "*(preview unavailable)*"
        if first_insert and Path(first_insert).exists():
            preview = Path(first_insert).read_text().splitlines()[0][:80] + "..."

        summary = textwrap.dedent(
            """
            ### Insert Profiling
            Extracted inserts and generated QC metrics. Download the archive for full outputs.
            """
        )
        archive = _zip_paths([Path(r["directory"]) for r in results], "profile_inserts")
        return summary, str(archive)
    except Exception as exc:  # pragma: no cover
        _LOGGER.exception("Profile inserts GUI failure")
        return f"⚠️ Error: {exc}", None


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
                    A guided graphical interface for primer design and sequencing analysis.
                    Use the tabs below, supply the required inputs, and download the generated results.
                    """
                )
            )

        with gr.Tab("Nextera XT"):  # --- Nextera ---
            gr.Markdown(
                """
                ### Illumina-Compatible Primer Design
                Provide the forward and reverse binding regions in 5'→3' orientation.
                """
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

        with gr.Tab("SLIM"):
            gr.Markdown(
                """
                ### Sequence-Ligation Independent Mutagenesis
                Paste the gene coding sequence, the plasmid context, and one mutation per line.
                """
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

        with gr.Tab("Gibson"):
            gr.Markdown(
                """
                ### Gibson Assembly Primer Design
                Use `+` to combine multiple mutations applied simultaneously.
                """
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

        with gr.Tab("Mutation Caller"):
            gr.Markdown(
                """
                ### Long-read Mutation Analysis
                Upload a FASTQ(.gz), the template FASTA, and the mutation_caller CSV configuration.
                """
            )
            mc_fastq = gr.File(label="FASTQ (.fastq.gz)", file_types=[".fastq", ".gz"], type="filepath")
            mc_template = gr.File(label="Template FASTA", file_types=[".fasta", ".fa"], type="filepath")
            mc_config = gr.File(label="Configuration CSV", file_types=[".csv"], type="filepath")
            mc_btn = gr.Button("Run mutation caller", variant="primary")
            mc_summary = gr.Markdown(label="Summary")
            mc_download = gr.File(label="Download results", file_count="single")
            mc_btn.click(
                fn=run_gui_mutation_caller,
                inputs=[mc_fastq, mc_template, mc_config],
                outputs=[mc_summary, mc_download],
            )

        with gr.Tab("UMI Hunter"):
            gr.Markdown(
                """
                ### UMI-Gene Pair Clustering
                Upload a FASTQ(.gz), template FASTA, and the UMI configuration CSV.
                """
            )
            umi_fastq = gr.File(label="FASTQ (.fastq.gz)", file_types=[".fastq", ".gz"], type="filepath")
            umi_template = gr.File(label="Template FASTA", file_types=[".fasta", ".fa"], type="filepath")
            umi_config = gr.File(label="UMI config CSV", file_types=[".csv"], type="filepath")
            umi_btn = gr.Button("Run UMI hunter", variant="primary")
            umi_summary = gr.Markdown(label="Summary")
            umi_download = gr.File(label="Download results", file_count="single")
            umi_btn.click(
                fn=run_gui_umi_hunter,
                inputs=[umi_fastq, umi_template, umi_config],
                outputs=[umi_summary, umi_download],
            )

        with gr.Tab("Profile Inserts"):
            gr.Markdown(
                """
                ### Insert Profiling
                Upload the probe CSV and one or more FASTQ(.gz) files containing reads.
                """
            )
            pi_csv = gr.File(label="Probe CSV", file_types=[".csv"], type="filepath")
            pi_fastq = gr.File(
                label="FASTQ files",
                file_types=[".fastq", ".gz"],
                file_count="multiple",
                type="filepath",
            )
            pi_btn = gr.Button("Profile inserts", variant="primary")
            pi_summary = gr.Markdown(label="Summary")
            pi_download = gr.File(label="Download results", file_count="single")
            pi_btn.click(
                fn=run_gui_profile_inserts,
                inputs=[pi_csv, pi_fastq],
                outputs=[pi_summary, pi_download],
            )

        with gr.Tab("EP Library Profile"):
            gr.Markdown(
                """
                ### Library Profiling Without UMIs
                Upload one or more FASTQ(.gz) files plus the region and plasmid references.
                """
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

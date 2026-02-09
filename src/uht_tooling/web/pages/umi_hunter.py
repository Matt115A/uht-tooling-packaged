"""UMI Hunter page — 11 inputs including sliders (most complex form)."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_input,
    apple_markdown,
    apple_number,
    apple_progress,
    apple_slider,
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_umi_hunter


async def render() -> None:
    with apple_card(
        "UMI Hunter",
        "Detect UMI barcodes, extract paired gene inserts, cluster reads by UMI "
        "identity, and emit consensus sequences with abundance tables.",
    ):
        # File state
        fastq_path: dict[str, Optional[str]] = {"value": None}
        template_path: dict[str, Optional[str]] = {"value": None}

        async def _save_upload(e, target: dict) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_umi_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            target["value"] = str(dest)

        async def _save_fastq(e):
            await _save_upload(e, fastq_path)

        async def _save_template(e):
            await _save_upload(e, template_path)

        # File uploads
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                apple_upload(
                    "FASTQ file (.fastq/.gz)",
                    extensions=[".fastq", ".gz"],
                    on_upload=_save_fastq,
                )
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Template FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_template,
                )

        template_text = apple_textarea(
            "Template FASTA (paste)",
            "Paste FASTA or raw DNA sequence here",
            rows=4,
        )

        def _normalize_fasta(text: str, header: str) -> str:
            cleaned = text.strip()
            if not cleaned:
                return ""
            if cleaned.lstrip().startswith(">"):
                return cleaned.rstrip() + "\n"
            seq = "".join(cleaned.split())
            lines = [seq[i : i + 80] for i in range(0, len(seq), 80)]
            return f">{header}\n" + "\n".join(lines) + "\n"

        def _write_fasta_from_text(text: str, header: str) -> Optional[str]:
            fasta = _normalize_fasta(text, header)
            if not fasta:
                return None
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_umi_"))
            dest = tmp / "pasted_template.fasta"
            dest.write_text(fasta)
            return str(dest)

        # UMI flanks
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                umi_start = apple_input(
                    "UMI upstream flank (5'\u21923')",
                    "e.g. ACACTCTTTCCCTACACGAC",
                )
            with ui.column().classes("flex-1"):
                umi_end = apple_input(
                    "UMI downstream flank (5'\u21923')",
                    "e.g. GACTGGAGTTCAGACGTGTG",
                )

        # Gene flanks
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                gene_start = apple_input(
                    "Gene upstream flank (5'\u21923')", "e.g. ATG..."
                )
            with ui.column().classes("flex-1"):
                gene_end = apple_input(
                    "Gene downstream flank (5'\u21923')", "e.g. TTA..."
                )

        # UMI length bounds
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                umi_min = apple_number("Minimum UMI length (nt)", 8, min_val=1)
            with ui.column().classes("flex-1"):
                umi_max = apple_number("Maximum UMI length (nt)", 14, min_val=1)

        # Sliders
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                umi_identity = apple_slider(
                    "UMI clustering identity",
                    min_val=0.5,
                    max_val=1.0,
                    value=0.9,
                    step=0.05,
                    unit="%",
                    display_multiplier=100.0,
                    precision=0,
                )
            with ui.column().classes("flex-1"):
                consensus_thresh = apple_slider(
                    "Consensus mutation threshold",
                    min_val=0.5,
                    max_val=1.0,
                    value=0.7,
                    step=0.05,
                    unit="%",
                    display_multiplier=100.0,
                    precision=0,
                )

        min_cluster = apple_slider(
            "Minimum reads per cluster",
            min_val=1,
            max_val=50,
            value=3,
            step=1,
            unit="reads",
            precision=0,
        )

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            template_value = template_path["value"]
            pasted_template = _write_fasta_from_text(
                template_text.value or "", "pasted_template"
            )
            if pasted_template:
                template_value = pasted_template

            summary, zip_path = await run_in_threadpool(
                run_gui_umi_hunter,
                fastq_path["value"],
                template_value,
                umi_start.value,
                umi_end.value,
                umi_min.value,
                umi_max.value,
                gene_start.value,
                gene_end.value,
                umi_identity.value,
                consensus_thresh.value,
                int(min_cluster.value),
            )

            progress.set_visibility(False)
            result_md.set_content(summary)

            if zip_path:
                download_row.clear()
                with download_row:
                    ui.button(
                        "Download Results",
                        on_click=lambda: ui.download(zip_path),
                    ).classes("apple-btn-secondary").props("unelevated no-caps")
                download_row.set_visibility(True)

        apple_button("Run UMI Hunter", on_click=on_run)

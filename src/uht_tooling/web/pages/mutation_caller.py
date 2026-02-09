"""Mutation Caller page — file uploads + number inputs."""

from __future__ import annotations

import shutil
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
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_mutation_caller


async def render() -> None:
    with apple_card(
        "Mutation Caller",
        "Extract coding regions bounded by user-defined flanks, align to the "
        "template, and report amino-acid substitutions.",
    ):
        # File state
        fastq_path: dict[str, Optional[str]] = {"value": None}
        template_path: dict[str, Optional[str]] = {"value": None}

        async def _save_upload(e, target: dict, suffix: str = "") -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_mc_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            target["value"] = str(dest)

        async def _save_fastq(e):
            await _save_upload(e, fastq_path)

        async def _save_template(e):
            await _save_upload(e, template_path)

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                fastq_upload = apple_upload(
                    "FASTQ file (.fastq/.gz)",
                    extensions=[".fastq", ".gz"],
                    on_upload=_save_fastq,
                )
            with ui.column().classes("flex-1"):
                template_upload = apple_upload(
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
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_mc_"))
            dest = tmp / "pasted_template.fasta"
            dest.write_text(fasta)
            return str(dest)

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                upstream = apple_input(
                    "Upstream flank (5'\u21923')", "e.g. ACTGTTAG"
                )
            with ui.column().classes("flex-1"):
                downstream = apple_input(
                    "Downstream flank (5'\u21923')", "e.g. CGAACCTA"
                )

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                min_len = apple_number("Minimum gene length (nt)", 900, min_val=1)
            with ui.column().classes("flex-1"):
                max_len = apple_number("Maximum gene length (nt)", 1200, min_val=1)

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
                run_gui_mutation_caller,
                fastq_path["value"],
                template_value,
                upstream.value,
                downstream.value,
                min_len.value,
                max_len.value,
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

        apple_button("Run Mutation Caller", on_click=on_run)

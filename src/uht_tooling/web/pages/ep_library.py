"""EP Library Profile page — multi-file upload with long-running workflow."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List, Optional

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_markdown,
    apple_progress,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_ep_library_profile


async def render() -> None:
    with apple_card(
        "EP Library Profile",
        "Estimate background and target mutation rates for enzyme evolution "
        "libraries without UMI barcodes. Upload one or more FASTQ files; "
        "each is profiled individually.",
    ):
        # File state
        fastq_paths: List[str] = []
        region_path: dict[str, Optional[str]] = {"value": None}
        plasmid_path: dict[str, Optional[str]] = {"value": None}

        async def _save_fastq(e) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ep_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            fastq_paths.append(str(dest))

        async def _save_single(e, target: dict) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ep_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            target["value"] = str(dest)

        apple_upload(
            "FASTQ file(s) (.fastq/.gz) \u2014 upload multiple to profile each individually",
            extensions=[".fastq", ".gz"],
            multiple=True,
            on_upload=_save_fastq,
        )

        async def _save_region(e):
            await _save_single(e, region_path)

        async def _save_plasmid(e):
            await _save_single(e, plasmid_path)

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Region-of-interest FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_region,
                )
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Plasmid FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_plasmid,
                )

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            summary, zip_path = await run_in_threadpool(
                run_gui_ep_library_profile,
                fastq_paths,
                region_path["value"],
                plasmid_path["value"],
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

        apple_button("Run Profiling", on_click=on_run)

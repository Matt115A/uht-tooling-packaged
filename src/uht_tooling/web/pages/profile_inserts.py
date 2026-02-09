"""Profile Inserts page — editable probe table + multi-file upload."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List, Optional

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_button_secondary,
    apple_card,
    apple_input,
    apple_markdown,
    apple_progress,
    apple_slider,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_profile_inserts


async def render() -> None:
    with apple_card(
        "Profile Inserts",
        "Characterise inserts demarcated by upstream/downstream probes, extract "
        "sequences, and produce QC plots plus summary tables.",
    ):
        # Editable probe pairs
        ui.label("Probe Pairs").classes("text-xs font-semibold").style(
            "color: var(--text-secondary); margin-bottom: 4px;"
        )

        probe_entries: List[dict] = []
        probe_container = ui.column().classes("w-full gap-2")

        def _create_probe_row(name: str = "") -> dict:
            """Add a row of ui.input fields and return the entry dict."""
            with probe_container:
                with ui.row().classes("w-full gap-4 items-end"):
                    name_input = apple_input("Name (optional)", value=name)
                    upstream_input = apple_input(
                        "Upstream (5'\u21923')", "ATCG..."
                    )
                    downstream_input = apple_input(
                        "Downstream (5'\u21923')", "ATCG..."
                    )
            entry = {
                "name": name_input,
                "upstream": upstream_input,
                "downstream": downstream_input,
            }
            probe_entries.append(entry)
            return entry

        _create_probe_row("probe_1")

        def add_probe_row() -> None:
            idx = len(probe_entries) + 1
            _create_probe_row(f"probe_{idx}")

        apple_button_secondary("+ Add Probe Row", on_click=add_probe_row)

        # Multi-file upload
        uploaded_paths: List[str] = []

        async def _save_fastq(e) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_pi_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            uploaded_paths.append(str(dest))

        apple_upload(
            "FASTQ files (.fastq/.gz)",
            extensions=[".fastq", ".gz"],
            multiple=True,
            on_upload=_save_fastq,
        )

        # Similarity slider
        min_ratio = apple_slider(
            "Minimum fuzzy-match ratio",
            min_val=50,
            max_val=100,
            value=80,
            step=1,
        )

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            import pandas as pd

            df = pd.DataFrame([
                {
                    "name": entry["name"].value,
                    "upstream": entry["upstream"].value,
                    "downstream": entry["downstream"].value,
                }
                for entry in probe_entries
            ])

            summary, zip_path = await run_in_threadpool(
                run_gui_profile_inserts,
                df,
                uploaded_paths,
                int(min_ratio.value),
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

        apple_button("Profile Inserts", on_click=on_run)

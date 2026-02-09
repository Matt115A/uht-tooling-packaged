"""Nextera XT primer design page."""

from __future__ import annotations

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_input,
    apple_markdown,
    apple_progress,
)
from uht_tooling.workflows.gui import run_gui_nextera


async def render() -> None:
    with apple_card(
        "Nextera XT Primer Design",
        "Generate Nextera XT-ready primers from forward/reverse binding regions.",
    ):
        fwd = apple_input("Forward primer (5'\u21923')", "e.g. ATCGATCG...")
        rev = apple_input("Reverse primer (5'\u21923')", "e.g. GCTAGCTA...")

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            summary, zip_path = await run_in_threadpool(
                run_gui_nextera, fwd.value, rev.value
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

        apple_button("Generate Primers", on_click=on_run)

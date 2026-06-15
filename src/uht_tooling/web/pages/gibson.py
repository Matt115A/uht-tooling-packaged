"""Gibson Assembly primer design page."""

from __future__ import annotations

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_markdown,
    apple_progress,
    apple_textarea,
)
from uht_tooling.workflows.gui import run_gui_design_gibson


async def render() -> None:
    with apple_card(
        "Gibson Assembly Primer Design",
        "Plan primer sets and assembly steps for Gibson mutagenesis, supporting "
        "multi-mutation constructs using the + syntax (e.g. A123G+T150A).",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow designs **Gibson assembly mutagenesis primers** and a matching assembly plan.

- You can specify **single mutations** or link multiple sub-mutations with `+` to define a combined construct, such as `A123G+T150A`.
- The tool breaks each construct into the PCR fragments needed for Gibson assembly and assigns the required primer pairs.
- It reports an explicit **assembly plan** with primer usage, fragment grouping, estimated Tm, and amplicon size.
- If requested mutations conflict within overlapping primer windows, the workflow can push you toward **sequential assembly steps** instead of a single reaction.

Outputs:

- `Gibson_primers.csv`: designed primers grouped by construct and submutation
- `Gibson_assembly_plan.csv`: PCR/assembly plan with primer assignments and amplicon sizes
                """
            ).classes("apple-markdown")
        gene = apple_textarea(
            "Gene sequence",
            "Paste the coding sequence for the gene of interest...",
            rows=4,
        )
        context = apple_textarea(
            "Plasmid context",
            "Paste the circular plasmid context sequence...",
            rows=4,
        )
        mutations = apple_textarea(
            "Mutations (one per line)",
            "e.g. A123G, A123G+T150A",
            rows=5,
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
                run_gui_design_gibson, gene.value, context.value, mutations.value
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

        apple_button("Design Gibson Primers", on_click=on_run)

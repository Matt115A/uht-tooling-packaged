"""SLIM primer design page."""

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
from uht_tooling.workflows.gui import run_gui_design_slim


async def render() -> None:
    with apple_card(
        "SLIM Primer Design",
        "Design paired short/long primers to introduce targeted mutations by "
        "Sequence-Ligation Independent Mutagenesis.",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow designs **SLIM mutagenesis primers** for targeted substitutions, deletions, insertions, and degenerate-codon libraries.

- For each mutation, the tool produces the standard **four-primer SLIM set**: long forward, short reverse, long reverse, and short forward.
- It validates each mutation against the supplied gene sequence before designing primers.
- Degenerate codons such as **`NNK`** or **`NNS`** are supported directly, so the same workflow can be used for site-saturation or focused library construction.
- Reverse primers are written with the correct **IUPAC reverse complements** for degenerate bases.

Mutation examples:

- `A123G` for a substitution
- `T241Del` for a deletion
- `T241TS` for an insertion
- `A123:NNK` for a degenerate library codon

Outputs:

- `SLIM_primers.csv`: four primers per mutation
                """
            ).classes("apple-markdown")
        gene = apple_textarea(
            "Gene sequence",
            "Paste the target gene coding sequence...",
            rows=4,
        )
        context = apple_textarea(
            "Plasmid context",
            "Paste the plasmid or genomic context...",
            rows=4,
        )
        mutations = apple_textarea(
            "Mutations (one per line)",
            "e.g. A123G, T241Del, R57:NNK",
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
                run_gui_design_slim, gene.value, context.value, mutations.value
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

        apple_button("Design SLIM Primers", on_click=on_run)

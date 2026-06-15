"""KLD (Inverse PCR) primer design page."""

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
from uht_tooling.workflows.gui import run_gui_design_kld


async def render() -> None:
    with apple_card(
        "KLD Primer Design",
        "Design inverse-PCR primers for KLD cloning. Forward primers carry "
        "the mutation at the 5\u2032 end; reverse primers bind upstream to "
        "re-amplify the full plasmid.",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow designs **KLD inverse-PCR primers** for introducing mutations into a plasmid with a two-primer strategy.

- For each mutation, the tool emits **one forward primer** and **one reverse primer**.
- The forward primer places the mutation at its **5' end** and extends downstream into the template-binding region.
- The reverse primer is designed immediately upstream so the pair can re-amplify the **entire plasmid**.
- Binding-region **Tm**, **GC%**, and **length** are reported for each primer, and the design aims to keep the two primer Tms closely matched.
- The same mutation syntax as SLIM is supported, including substitutions, deletions, insertions, and degenerate library codons.

Compared with SLIM, this is a simpler **two-primer** workflow and is typically better suited to **single-mutation** jobs.

Outputs:

- `KLD_primers.csv`: forward/reverse primer pairs plus binding metrics and notes
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
                run_gui_design_kld, gene.value, context.value, mutations.value
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

        apple_button("Design KLD Primers", on_click=on_run)

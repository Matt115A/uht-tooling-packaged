"""Gene oligo design page."""

from __future__ import annotations

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_markdown,
    apple_number,
    apple_progress,
    apple_textarea,
)
from uht_tooling.workflows.gui import (
    run_gui_design_gene_oligos,
    run_gui_length_optimize_gene_oligos,
)


async def render() -> None:
    with apple_card(
        "Cheap Genes",
        "Design overlap-extension PCR oligos for one or more IVTT-ready constructs from DNA or protein FASTA input.",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow designs **overlap-extension PCR oligos** for building full IVTT-ready gene constructs from shorter ordered DNA pieces.

- Each target is expanded into a full construct with the configured **5' constant region**, **3' constant region**, and any selected **N-terminal or C-terminal tag**.
- The construct is then tiled into **gene-specific overlap oligos** that can be assembled by overlap-extension PCR.
- The workflow also emits reusable **order-once constant oligos** at the construct edges, so you do not need to reorder the promoter/terminator-bearing pieces for every target.
- Separate **external primers** are produced for downstream full-length amplification.
- For **DNA inputs**, the tool removes a terminal stop codon when needed and adds a start codon if the construct requires one.
- For **protein inputs**, the sequence is first **codon-optimized for the selected host**, then converted into DNA before oligo tiling.
- The design engine chooses the **longest feasible oligo size at or below your requested maximum** so it can reduce primer count without violating overlap constraints.
- The **Length Optimize** action scans candidate oligo-length caps and ranks alternative assembly strategies, which is useful when you want to minimize oligo count or compare cost/complexity tradeoffs.

Outputs:

- `gene_oligos.csv`: all constant and gene-specific OE oligos with overlap metrics
- `assembly_report.csv`: per-target construct summary and warnings
- `ordered_oligos.fasta`: deduplicated ordering FASTA
- `ordering_plate.csv`: ordering-friendly plate layout
- `external_primers.csv`: reusable outer primers for full-length amplification
                """
            ).classes("apple-markdown")
        sequence = apple_textarea(
            "DNA or protein sequence(s)",
            "Paste multi-FASTA or a single raw sequence here...",
            rows=8,
        )
        with ui.row().classes("w-full gap-4"):
            input_type = ui.select(
                options={"auto": "Auto-detect", "dna": "DNA", "protein": "Protein"},
                value="auto",
                label="Input type",
            ).classes("apple-field flex-1")
            input_type.props("outlined dense")
            tag_mode = ui.select(
                options={
                    "none": "No tag",
                    "n_his": "N-terminal His",
                    "c_his": "C-terminal His",
                    "n_his_flag": "N-terminal His + FLAG",
                    "c_his_flag": "C-terminal His + FLAG",
                },
                value="none",
                label="Tag mode",
            ).classes("apple-field flex-1")
            tag_mode.props("outlined dense")

        target_oligo_length = ui.number(
            label="Max gene-specific oligo length (nt)",
            value=40,
            min=20,
            step=1,
            format="%.0f",
        ).classes("apple-field w-full")
        target_oligo_length.props("outlined dense")

        with ui.expansion("Extra settings", icon="settings").classes("w-full"):
            max_order_once_length = apple_number(
                "Max order-once oligo length (nt)",
                value=80,
                min_val=20,
                step=1,
            )
            target_host = ui.select(
                options={
                    "e_coli": "E. coli",
                    "b_subtilis": "B. subtilis",
                    "s_cerevisiae": "S. cerevisiae",
                    "p_pastoris": "P. pastoris",
                    "h_sapiens": "H. sapiens",
                    "m_musculus": "Mouse",
                    "c_elegans": "C. elegans",
                    "d_melanogaster": "D. melanogaster",
                    "a_thaliana": "A. thaliana",
                },
                value="e_coli",
                label="Target host for protein codon optimization",
            ).classes("apple-field w-full")
            target_host.props("outlined dense")

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            summary, zip_path = await run_in_threadpool(
                run_gui_design_gene_oligos,
                sequence.value,
                input_type.value,
                target_oligo_length.value,
                max_order_once_length.value,
                target_host.value,
                tag_mode.value,
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

        async def on_optimize() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            summary, zip_path = await run_in_threadpool(
                run_gui_length_optimize_gene_oligos,
                sequence.value,
                input_type.value,
                target_oligo_length.value,
                max_order_once_length.value,
                target_host.value,
                tag_mode.value,
            )

            progress.set_visibility(False)
            result_md.set_content(summary)

            if zip_path:
                download_row.clear()
                with download_row:
                    ui.button(
                        "Download Optimization Results",
                        on_click=lambda: ui.download(zip_path),
                    ).classes("apple-btn-secondary").props("unelevated no-caps")
                download_row.set_visibility(True)

        with ui.row().classes("w-full gap-3"):
            apple_button("Design Cheap Genes", on_click=on_run)
            ui.button("Length Optimize", on_click=on_optimize).classes("apple-btn-secondary").props("unelevated no-caps")

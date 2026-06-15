"""Synthetic gene pool ordering page."""

from __future__ import annotations

from pathlib import Path

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_markdown,
    apple_progress,
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_design_synthetic_gene_pool


async def render() -> None:
    with apple_card(
        "Synthetic Genes",
        "Create a pooled synthetic-gene ordering list for E. coli T7 cell-free protein synthesis, with shared lift-out primers that add the final 5' and 3' sequence features.",
    ):
        ui.image("/uht-static/animations/synthetic_gene_pool_workflow.gif").style(
            "width: 100%; border-radius: 12px; border: 1px solid var(--border);"
            " box-shadow: var(--shadow-md); margin-bottom: 8px; display: block;"
        )
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow designs a pooled ordering set for **E. coli T7 cell-free protein synthesis**.

- Each input target becomes **one pooled ordering oligo**.
- That pooled oligo carries the **shared anneal handles** plus the target coding region, so the ordered pool oligos stay shorter than a full expression cassette.
- The tool also emits **one reusable primer pair**: `POOL_CONST_F` and `POOL_CONST_R`.
- During lift-out PCR, those common primers add the final **5' constant region**, **3' constant region**, and any selected **N-terminal or C-terminal tag**.
- Optional **gene-specific pullout primers** can also be emitted. In that mode, each design gets a deterministic balanced index built into the ordered DNA just upstream of the terminator, and a matching reverse primer binds that unique region to selectively recover the design after whole-pool amplification.
- Intended pullout workflow:
  1. Order the pooled oligos, plus `POOL_CONST_F` and `POOL_CONST_R`.
  2. Use `POOL_CONST_F` + `POOL_CONST_R` to amplify the entire pool into a CFPS-ready library.
  3. Use `POOL_CONST_F` + the matching `*_PULLOUT_R` primer to selectively recover one design from that amplified pool.
  4. The recovered product already contains promoter, gene, unique index, and terminator, so it can be used directly for CFPS.
- For **protein inputs**, sequences are automatically **codon-optimized for E. coli** before the pooled oligos are written.
- Protein residue **`X`** is translated to **`NNN`** in the pooled oligo output so you can preserve a randomized codon at that position.
- The common primers must remain **below 100 bp**. If the configured constant regions plus tag payload are too long, the workflow stops with an error instead of producing an unusable ordering design.

Outputs:

- `synthetic_gene_pool.csv`: one pooled oligo per target
- `synthetic_gene_pool_primers.csv`: the shared lift-out primer pair to order once
- `synthetic_gene_pool_ordering.tsv`: combined ordering sheet for pool oligos and common primers
- `synthetic_gene_pool_instructions.txt`: brief wet-lab usage notes
                """
            ).classes("apple-markdown")
        ui.label(
            "For protein inputs, use `X` to denote a random amino acid. The pooled oligo output will use `NNN` at that codon position."
        ).classes("apple-card-subtitle w-full")
        sequence = apple_textarea(
            "DNA or protein sequence(s)",
            "Paste multi-FASTA or a single raw sequence here. Protein sequences may include X for random amino acids.",
            rows=8,
        )

        async def _append_uploaded_fasta(e) -> None:
            content = (await e.file.read()).decode("utf-8", errors="replace").strip()
            if not content:
                ui.notify(f"{e.file.name} was empty.", color="warning")
                return

            if not content.lstrip().startswith(">"):
                header = Path(e.file.name).stem or "uploaded_sequence"
                compact = "".join(content.split())
                lines = [compact[i : i + 80] for i in range(0, len(compact), 80)]
                content = f">{header}\n" + "\n".join(lines)

            existing = (sequence.value or "").strip()
            sequence.value = f"{existing}\n{content}\n" if existing else f"{content}\n"
            ui.notify(f"Added {e.file.name}", color="positive")

        apple_upload(
            "Upload FASTA file(s) (.fasta/.fa) — drag and drop one or more files to append them here",
            extensions=[".fasta", ".fa", ".faa"],
            multiple=True,
            on_upload=_append_uploaded_fasta,
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
        include_pullout_primers = ui.checkbox(
            "Add deterministic per-gene pullout primers for selective PCR recovery"
        ).classes("w-full")
        ui.label(
            "Pullout mode: first amplify the full pool with POOL_CONST_F/POOL_CONST_R, then use POOL_CONST_F plus a matching *_PULLOUT_R primer to recover one CFPS-ready design."
        ).classes("apple-card-subtitle w-full")

        ui.label(
            "Protein inputs are automatically codon-optimized for E. coli because this tool is specialized for T7-driven cell-free expression."
        ).classes("apple-card-subtitle w-full")

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            summary, zip_path = await run_in_threadpool(
                run_gui_design_synthetic_gene_pool,
                sequence.value,
                input_type.value,
                tag_mode.value,
                include_pullout_primers.value,
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

        apple_button("Design Synthetic Gene Pool", on_click=on_run)

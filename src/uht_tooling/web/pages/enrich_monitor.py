"""Enrichment Monitor page — variant enrichment quantification from Nanopore reads."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List, Optional

from starlette.concurrency import run_in_threadpool
from nicegui import ui

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_input,
    apple_markdown,
    apple_progress,
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_enrich_monitor


async def render() -> None:
    with apple_card(
        "Enrichment Monitor",
        "Quantify variant enrichment from Nanopore long-read sequencing data — "
        "map reads to a multi-sequence reference and compute per-sample proportions, "
        "Wilson confidence intervals, and Fisher exact enrichment odds ratios.",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow quantifies **how much a target sequence is enriched** relative to a background
(non-binding) control sequence after one or more rounds of affinity selection or partitioning.

**Typical use-cases:**
- Yeast surface display (YSD) enrichment monitoring between selection rounds
- Phage display panning: track target phage vs. helper phage ratio
- Any experiment where you have a mixed pool of two or more sequence variants and
  want to know which is enriched in the output vs. the input

**How it works:**
1. Reads from each FASTQ are aligned to the reference FASTA (containing all variants) using minimap2.
2. Primary alignments are counted per reference record.
3. The target fraction and 95 % Wilson CI are computed for each sample.
4. A Fisher exact odds ratio and Woolf 95 % CI are computed for each selected sample
   relative to the designated unselected baseline.
5. Results are exported as TSV tables and a two-panel summary figure.

**Reference FASTA format:**
Two or more records, one per variant. Example:
```
>CBM1
ATGAAGAAA...
>Helix_alone
ATGGCAAAA...
```

The positive control (target) and negative control (background) sequences can be fully different
genes, point mutants, or deletion variants — anything distinguishable by Nanopore alignment.
                """
            ).classes("apple-markdown")

        # ── State ─────────────────────────────────────────────────────────────
        fastq_paths: List[str] = []
        ref_path: dict[str, Optional[str]] = {"value": None}

        async def _save_fastq(e) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_enrich_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            dest.write_bytes(data)
            fastq_paths.append(str(dest))

        async def _save_ref(e) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_enrich_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            dest.write_bytes(data)
            ref_path["value"] = str(dest)

        # ── Reference FASTA ────────────────────────────────────────────────────
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Reference FASTA (.fasta/.fa) — upload",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_ref,
                )
                ui.label(
                    "Multi-sequence FASTA containing all variant sequences to count against "
                    "(target + background, and any additional variants)."
                ).classes("apple-card-subtitle w-full")
            with ui.column().classes("flex-1"):
                ref_text = apple_textarea(
                    "Reference FASTA — paste",
                    "Paste multi-sequence FASTA here (≥2 records)",
                    rows=5,
                )

        # ── FASTQ files ────────────────────────────────────────────────────────
        apple_upload(
            "FASTQ file(s) (.fastq/.fastq.gz) — upload all conditions including the unselected baseline",
            extensions=[".fastq", ".gz"],
            multiple=True,
            on_upload=_save_fastq,
        )
        ui.label(
            "Upload all samples: both the unselected input pool and every selected output condition. "
            "Each file is treated as one independent sample; the filename stem (without .fastq.gz) "
            "becomes the sample name."
        ).classes("apple-card-subtitle w-full")

        # ── Parameter inputs ───────────────────────────────────────────────────
        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                target_input = apple_input(
                    "Target sequence ID",
                    "e.g. CBM1  (defaults to first record)",
                )
            with ui.column().classes("flex-1"):
                background_input = apple_input(
                    "Background sequence ID",
                    "e.g. Helix_alone  (defaults to second record)",
                )

        baseline_input = apple_input(
            "Baseline sample name (unselected input pool)",
            "e.g. control  — must match the FASTQ file stem exactly",
        )

        with ui.row().classes("w-full gap-4 items-center"):
            n_bootstrap_input = ui.number(
                label="Bootstrap iterations (Baret η / Zinchenko η' CI)",
                value=2000,
                min=100,
                max=50000,
                step=500,
            ).classes("flex-none w-72")
            ui.label(
                "Number of binomial resampling draws used to estimate 95 % CIs for the "
                "Baret odds-ratio and Zinchenko relative-risk enrichment scores."
            ).classes("apple-card-subtitle flex-1")

        # ── Run ────────────────────────────────────────────────────────────────
        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            ref_fasta_content = ""
            if ref_path["value"]:
                ref_fasta_content = Path(ref_path["value"]).read_text()
            elif ref_text.value:
                ref_fasta_content = ref_text.value

            summary, zip_path = await run_in_threadpool(
                run_gui_enrich_monitor,
                fastq_paths,
                ref_fasta_content,
                target_input.value or "",
                background_input.value or "",
                baseline_input.value or "",
                int(n_bootstrap_input.value or 2000),
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

        apple_button("Run Enrichment Analysis", on_click=on_run)

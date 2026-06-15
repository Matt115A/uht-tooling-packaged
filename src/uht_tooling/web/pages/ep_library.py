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
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_ep_library_profile


async def render() -> None:
    with apple_card(
        "EP Library Profile",
        "Quantify cloning accuracy and mutation rates in a gene of interest, "
        "identify commonly mutated sites in a plasmid pool and their relative "
        "abundances, and compare variant enrichments between DNA pools.",
    ):
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow profiles **error-prone or heterogeneous plasmid pools** from long-read sequencing data without UMIs.

- It aligns reads to both the **region of interest (ROI)** and the full **plasmid**.
- The plasmid is treated as **circular**, so the ROI can be forward, reverse-complemented, or split across the plasmid origin and still be found.
- Mismatch rates inside the ROI are compared against plasmid-wide background outside the ROI to estimate a **net mutation rate**.
- The profiler reports a single **lambda** value, interpreted as mutations per gene copy, and uses simulation to estimate the expected **amino-acid mutation burden**.
- Each sample gets both technical outputs and a more readable summary package.

Outputs include:

- per-sample summary panels
- detailed mismatch and amino-acid substitution tables
- sample reports and mutation-spectrum plots
- top-level master summary when multiple FASTQs are processed
                """
            ).classes("apple-markdown")
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
        ui.label(
            "Raw outputs of a nanopore sequencing run, for example from a whole-plasmid sequencing run."
        ).classes("apple-card-subtitle w-full")

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
                ui.label(
                    "FASTA of the gene of interest, from start codon to stop codon."
                ).classes("apple-card-subtitle w-full")
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Plasmid FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_plasmid,
                )
                ui.label(
                    "FASTA of the plasmid containing the gene of interest. The plasmid is treated as circular, so the gene may be forward, reverse-complemented, or split across the start and end of the entry, but it must still match uniquely."
                ).classes("apple-card-subtitle w-full")

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                region_text = apple_textarea(
                    "Region-of-interest FASTA (paste)",
                    "Paste FASTA or raw DNA sequence here",
                    rows=4,
                )
            with ui.column().classes("flex-1"):
                plasmid_text = apple_textarea(
                    "Plasmid FASTA (paste)",
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

        def _write_fasta_from_text(
            text: str, header: str, filename: str
        ) -> Optional[str]:
            fasta = _normalize_fasta(text, header)
            if not fasta:
                return None
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ep_"))
            dest = tmp / filename
            dest.write_text(fasta)
            return str(dest)

        progress = apple_progress()
        result_md = apple_markdown()
        download_row = ui.element("div")
        download_row.set_visibility(False)

        async def on_run() -> None:
            progress.set_visibility(True)
            result_md.set_content("")
            download_row.set_visibility(False)

            region_value = region_path["value"]
            plasmid_value = plasmid_path["value"]
            pasted_region = _write_fasta_from_text(
                region_text.value or "", "pasted_region", "pasted_region.fasta"
            )
            pasted_plasmid = _write_fasta_from_text(
                plasmid_text.value or "", "pasted_plasmid", "pasted_plasmid.fasta"
            )
            if pasted_region:
                region_value = pasted_region
            if pasted_plasmid:
                plasmid_value = pasted_plasmid

            summary, zip_path = await run_in_threadpool(
                run_gui_ep_library_profile,
                fastq_paths,
                region_value,
                plasmid_value,
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

"""SSM Profiler page."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List, Optional

from nicegui import ui
from starlette.concurrency import run_in_threadpool

from uht_tooling.web.components import (
    apple_button,
    apple_card,
    apple_markdown,
    apple_progress,
    apple_textarea,
    apple_upload,
)
from uht_tooling.workflows.gui import run_gui_ssm_profiler


async def render() -> None:
    with apple_card(
        "SSM Profiler",
        "Profile site-saturation libraries at user-defined codons, quantify observed amino-acid "
        "distributions, compare them to degenerate-codon expectations, and assess off-target "
        "mutation elsewhere in the gene of interest.",
    ):
        ui.image("/uht-static/animations/ssm_profiler.gif").style(
            "width: 100%; border-radius: 12px; border: 1px solid var(--border);"
            " box-shadow: var(--shadow-md); margin-bottom: 8px; display: block;"
        )
        with ui.expansion("What This Tool Does", icon="info").classes("w-full").props("default-opened"):
            ui.markdown(
                """
This workflow profiles **site-saturation mutagenesis (SSM) libraries** at user-specified codons from long-read sequencing data.

- You define the **ROI coding sequence**, the full **plasmid**, and the amino-acid positions that were intentionally diversified.
- Optionally, you can provide the **degenerate codon scheme** used at each site, such as `45:NNK`.
- The tool reports **observed amino-acid distributions** at each target site and, when schemes are supplied, compares them with the **expected distributions** from the degenerate code.
- Target-site mutational load is computed using only the codons you specify.
- Mismatches elsewhere in the ROI are treated as **off-target signal** and compared against plasmid background outside the ROI.

Outputs focus on:

- per-site amino-acid distributions
- target-only mutational-load summaries
- observed-vs-expected comparisons for supplied degenerate schemes
- off-target mismatch rates elsewhere in the ROI
                """
            ).classes("apple-markdown")

        fastq_paths: List[str] = []
        region_path: dict[str, Optional[str]] = {"value": None}
        plasmid_path: dict[str, Optional[str]] = {"value": None}

        async def _save_fastq(e) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ssm_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            fastq_paths.append(str(dest))

        async def _save_single(e, target: dict[str, Optional[str]]) -> None:
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ssm_"))
            dest = tmp / e.file.name
            data = await e.file.read()
            with open(dest, "wb") as f:
                f.write(data)
            target["value"] = str(dest)

        async def _save_region(e) -> None:
            await _save_single(e, region_path)

        async def _save_plasmid(e) -> None:
            await _save_single(e, plasmid_path)

        apple_upload(
            "FASTQ file(s) (.fastq/.gz) — upload multiple to profile each individually",
            extensions=[".fastq", ".gz"],
            multiple=True,
            on_upload=_save_fastq,
        )
        ui.label(
            "Raw long-read sequencing files, for example from a whole-plasmid nanopore run."
        ).classes("apple-card-subtitle w-full")

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                apple_upload(
                    "ROI CDS FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_region,
                )
                ui.label(
                    "FASTA of the coding sequence being diversified, from start codon to stop codon. "
                    "Target-site numbering is relative to this CDS."
                ).classes("apple-card-subtitle w-full")
            with ui.column().classes("flex-1"):
                apple_upload(
                    "Plasmid FASTA (.fasta/.fa)",
                    extensions=[".fasta", ".fa"],
                    on_upload=_save_plasmid,
                )
                ui.label(
                    "FASTA of the full plasmid containing the ROI. The plasmid alignment is used to estimate "
                    "background mismatches outside the ROI."
                ).classes("apple-card-subtitle w-full")

        with ui.row().classes("w-full gap-4"):
            with ui.column().classes("flex-1"):
                region_text = apple_textarea(
                    "ROI CDS FASTA (paste)",
                    "Paste FASTA or raw coding DNA sequence here",
                    rows=4,
                )
            with ui.column().classes("flex-1"):
                plasmid_text = apple_textarea(
                    "Plasmid FASTA (paste)",
                    "Paste FASTA or raw DNA sequence here",
                    rows=4,
                )

        target_sites = apple_textarea(
            "Target AA sites",
            "Comma-separated amino-acid positions, e.g. 45,46,47",
            rows=3,
        )
        ui.label(
            "Use amino-acid numbering relative to the ROI CDS above. Enter only the intentionally diversified codons."
        ).classes("apple-card-subtitle w-full")
        site_schemes = apple_textarea(
            "Optional site schemes",
            "One per line, e.g. 45:NNK",
            rows=4,
        )
        ui.label(
            "Optional. Map each target site to the degenerate codon used during library construction, "
            "for example `45:NNK`, `46:NNK`, `47:NNW`."
        ).classes("apple-card-subtitle w-full")

        ui.markdown(
            """
            **Interpretation notes**

            - The target-site mutational load is computed only from the codons you specify.
            - Off-target signal means mismatches elsewhere in the ROI, outside those target codons.
            - Off-target ROI mismatch rates are compared against plasmid-wide background outside the ROI.
            """
        ).classes("apple-card-subtitle w-full")

        def _normalize_fasta(text: str, header: str) -> str:
            cleaned = text.strip()
            if not cleaned:
                return ""
            if cleaned.lstrip().startswith(">"):
                return cleaned.rstrip() + "\n"
            seq = "".join(cleaned.split())
            lines = [seq[i : i + 80] for i in range(0, len(seq), 80)]
            return f">{header}\n" + "\n".join(lines) + "\n"

        def _write_fasta_from_text(text: str, header: str, filename: str) -> Optional[str]:
            fasta = _normalize_fasta(text, header)
            if not fasta:
                return None
            tmp = Path(tempfile.mkdtemp(prefix="uht_gui_ssm_"))
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
                region_text.value or "", "pasted_roi_cds", "pasted_roi_cds.fasta"
            )
            pasted_plasmid = _write_fasta_from_text(
                plasmid_text.value or "", "pasted_plasmid", "pasted_plasmid.fasta"
            )
            if pasted_region:
                region_value = pasted_region
            if pasted_plasmid:
                plasmid_value = pasted_plasmid

            summary, zip_path = await run_in_threadpool(
                run_gui_ssm_profiler,
                fastq_paths,
                region_value,
                plasmid_value,
                target_sites.value or "",
                site_schemes.value or "",
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

        apple_button("Run SSM Profiling", on_click=on_run)

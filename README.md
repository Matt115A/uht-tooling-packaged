# uht-tooling

Automation helpers for ultra-high-throughput molecular biology workflows. The package ships both a CLI and an optional GUI that wrap the same workflow code paths.

---

## Installation

### Quick install (recommended, easiest file maintainance)
```bash
pip install "uht-tooling[gui]"

```

This installs the core workflows plus the optional GUI dependency (Gradio). Omit the `[gui]` extras if you only need the CLI:

```bash
pip install uht-tooling
```

You will need a functioning version of mafft - you should install this separately and it should be accessible from your environment.

### Development install
```bash
git clone https://github.com/Matt115A/uht-tooling-packaged.git
cd uht-tooling-packaged
python -m pip install -e ".[gui,dev]"
```

The editable install exposes the latest sources, while the `dev` extras add linting and test tooling.

---

## Directory layout

- Reference inputs can be found anywhere (you specify in the cli), but we recommend using `data/<workflow>/`.
- Outputs (CSV, FASTA, plots, logs) are written to `results/<workflow>/`.
- All workflows log to `results/<workflow>/run.log` for reproducibility and debugging.

---

## Command-line interface

The CLI is exposed as the `uht-tooling` executable. List the available commands:

```bash
uht-tooling --help
```

Each command mirrors a workflow module. Common entry points:

| Command | Purpose |
| --- | --- |
| `uht-tooling nextera-primers` | Generate Nextera XT primer pairs from a binding-region CSV. |
| `uht-tooling design-slim` | Design SLIM mutagenesis primers from FASTA/CSV inputs. |
| `uht-tooling design-kld` | Design KLD (inverse PCR) mutagenesis primers. |
| `uht-tooling design-gibson` | Produce Gibson mutagenesis primers and assembly plans. |
| `uht-tooling mutation-caller` | Summarise amino-acid substitutions from long-read FASTQ files. |
| `uht-tooling umi-hunter` | Cluster UMIs and call consensus genes. |
| `uht-tooling ep-library-profile` | Measure mutation rates in plasmid libraries without UMIs. |
| `uht-tooling profile-inserts` | Extract and analyse inserts defined by flanking probe pairs. |

Each command provides detailed help, including option descriptions and expected file formats:

```bash
uht-tooling mutation-caller --help
```

You can pass multiple FASTQ paths using repeated `--fastq` options or glob patterns. Optional `--log-path` flags redirect logs if you prefer a location outside the default results directory.

---

## Workflow reference

### Nextera XT primer design

1. Prepare `data/nextera_designer/nextera_designer.csv` with a `binding_region` column. Row 1 should contain the forward region, row 2 the reverse region, both in 5'→3' orientation.
2. Optional: supply a YAML overrides file for index lists/prefixes via `--config`.
3. Run:
   ```bash
   uht-tooling nextera-primers \
     --binding-csv data/nextera_designer/nextera_designer.csv \
     --output-csv results/nextera_designer/nextera_xt_primers.csv
   ```
4. Primer CSVs will be written to `results/nextera_designer/`, accompanied by a log file.

The helper is preloaded with twelve i5 and twelve i7 indices, enabling up to 144 unique amplicons.

#### Wet-lab workflow notes

- Perform the initial amplification with an i5/i7 primer pair and monitor a small aliquot by qPCR. Cap thermocycling early so you only generate ~10% of the theoretical yield—this minimizes amplification bias.
- Purify the product with SPRIselect beads at approximately a 0.65:1 bead:DNA volume ratio to remove residual primers and short fragments.
- Confirm primer removal and quantify DNA using electrophoresis (e.g., BioAnalyzer DNA chip) before moving to the flow cell.

### SLIM primer design

- Inputs:
  - `data/design_slim/slim_template_gene.fasta`
  - `data/design_slim/slim_context.fasta`
  - `data/design_slim/slim_target_mutations.csv` (single `mutations` column)
- Run:
  ```bash
  uht-tooling design-slim \
    --gene-fasta data/design_slim/slim_template_gene.fasta \
    --context-fasta data/design_slim/slim_context.fasta \
    --mutations-csv data/design_slim/slim_target_mutations.csv \
    --output-dir results/design_slim/
  ```
- Output: `results/design_slim/SLIM_primers.csv` plus logs.

Mutation nomenclature examples:
- `A123G` (substitution)
- `T241Del` (deletion)
- `T241TS` (insert Ser after Thr241)
- `L46GP` (replace Leu46 with Gly-Pro)
- `A123:NNK` (library mutation with degenerate codon)

#### Library mutations with degenerate codons

For saturation mutagenesis and library generation, SLIM supports degenerate (IUPAC ambiguity) codons using the format `<WT_AA><position>:<codon>`. The codon must be exactly 3 characters using valid IUPAC nucleotide codes:

| Code | Bases | Mnemonic |
|------|-------|----------|
| A, C, G, T | Single base | Standard |
| R | A, G | puRine |
| Y | C, T | pYrimidine |
| S | G, C | Strong |
| W | A, T | Weak |
| K | G, T | Keto |
| M | A, C | aMino |
| B | C, G, T | not A |
| D | A, G, T | not C |
| H | A, C, T | not G |
| V | A, C, G | not T |
| N | A, C, G, T | aNy |

Common degenerate codon schemes for library construction:

| Scheme | Codons | Amino acids | Stop codons | Notes |
|--------|--------|-------------|-------------|-------|
| NNK | 32 | 20 | 1 (TAG) | Reduced stop codon frequency |
| NNS | 32 | 20 | 1 (TAG) | Equivalent to NNK |
| NNN | 64 | 20 | 3 | All codons, higher stop frequency |
| NDT | 12 | 12 | 0 | F, L, I, V, Y, H, N, D, C, R, S, G only |

Example CSV with mixed mutation types:
```csv
mutations
A123G
T50:NNK
S100:NNS
T241Del
```

The workflow validates that the wild-type amino acid matches the template sequence and logs library coverage information (number of possible codons and amino acids) for each degenerate mutation. Primers are generated with the degenerate bases embedded; reverse primers contain the correct IUPAC reverse complements (e.g., K↔M, R↔Y, S↔S).

#### Experimental blueprint

- Hands-on time is approximately three hours (excluding protein purification), with mutant protein obtainable in roughly three days.
- Conduct two PCRs per mutant set: (A) long forward with short reverse and (B) long reverse with short forward.
- Combine 10 µL from each PCR with 10 µL H-buffer (150 mM Tris pH 8, 400 mM NaCl, 60 mM EDTA) for a 30 µL annealing reaction: 99 °C for 3 min, then two cycles of 65 °C for 5 min followed by 30 °C for 15 min, hold at 4 °C.
- Transform directly into NEB 5-alpha or BL21 (DE3) cells without additional cleanup. The protocol has been validated for simultaneous introduction of dozens of mutations.

### KLD primer design

KLD (Kinase-Ligation-DpnI) is an alternative mutagenesis method using inverse PCR to amplify the entire plasmid with mutations incorporated at the primer junction.

- Inputs: Same as SLIM design
  - `data/design_kld/kld_template_gene.fasta`
  - `data/design_kld/kld_context.fasta`
  - `data/design_kld/kld_target_mutations.csv` (single `mutations` column)
- Run:
  ```bash
  uht-tooling design-kld \
    --gene-fasta data/design_kld/kld_template_gene.fasta \
    --context-fasta data/design_kld/kld_context.fasta \
    --mutations-csv data/design_kld/kld_target_mutations.csv \
    --output-dir results/design_kld/
  ```
- Output: `results/design_kld/KLD_primers.csv` plus logs.

Mutation nomenclature: Same as SLIM (substitution, deletion, insertion, indel, library).

#### KLD vs SLIM

| Method | Primers | Mechanism | Best for |
|--------|---------|-----------|----------|
| SLIM | 4 per mutation | Overlap assembly | Multiple simultaneous mutations |
| KLD | 2 per mutation | Inverse PCR + ligation | Single mutations, simpler workflow |

#### KLD primer design rules

- Forward primer: Mutation codon at 5' end + downstream template-binding region
- Reverse primer: Reverse complement of upstream region, 5' end adjacent to forward
- Tm calculated on template-binding regions only (50-65°C target)
- Tm difference between primers kept within 5°C
- GC content 40-60%
- Binding region 18-24 bp

#### Experimental workflow

1. PCR amplify entire plasmid with KLD primer pair
2. DpnI digest to remove methylated template
3. T4 PNK phosphorylation of 5' ends
4. T4 DNA ligase to circularize
5. Transform into competent cells

NEB sells a KLD Enzyme Mix (M0554) that combines these steps.

### Gibson assembly primers

- Inputs mirror the SLIM workflow but use `data/design_gibson/`.
- Link sub-mutations with `+` to specify multi-mutation assemblies (e.g., `A123G+T150A`).
- Run:
  ```bash
  uht-tooling design-gibson \
    --gene-fasta data/design_gibson/gibson_template_gene.fasta \
    --context-fasta data/design_gibson/gibson_context.fasta \
    --mutations-csv data/design_gibson/gibson_target_mutations.csv \
    --output-dir results/design_gibson/
  ```
- Outputs include primer sets and an assembly-plan CSV.

If mutations fall within overlapping primer windows, design sequential reactions. 

### Mutation caller (no UMIs)

1. Supply:
   - `data/mutation_caller/mutation_caller_template.fasta`
   - `data/mutation_caller/mutation_caller.csv` with `gene_flanks` and `gene_min_max` columns (two rows each).
   - One or more FASTQ files via `--fastq`.
2. Run:
   ```bash
   uht-tooling mutation-caller \
     --template-fasta data/mutation_caller/mutation_caller_template.fasta \
     --flanks-csv data/mutation_caller/mutation_caller.csv \
     --fastq data/mutation_caller/*.fastq.gz \
     --output-dir results/mutation_caller/ \
     --threshold 10
   ```
3. Outputs: per-sample subdirectories with substitution summaries, co-occurrence matrices, and logs. Co-occurence matrices are experimental and are not yet to be relied on.

### UMI Hunter

- Inputs: `data/umi_hunter/template.fasta`, `data/umi_hunter/umi_hunter.csv`, and FASTQ reads.
- Command:
  ```bash
  uht-tooling umi-hunter \
    --template-fasta data/umi_hunter/template.fasta \
    --config-csv data/umi_hunter/umi_hunter.csv \
    --fastq data/umi_hunter/*.fastq.gz \
    --output-dir results/umi_hunter/
  ```
- Tunable parameters include `--umi-identity-threshold`, `--consensus-mutation-threshold`, and `--min-cluster-size`.
- `--umi-identity-threshold` (0–1) controls how similar two UMIs must be to fall into the same cluster.
- `--consensus-mutation-threshold` (0–1) is the fraction of reads within a cluster that must agree on a base before it is written into the consensus sequence.
- `--min-cluster-size` sets the minimum number of reads required in a cluster before a consensus is generated (smaller clusters remain listed in the raw UMI CSV but no consensus FASTA is produced).

Please be aware, this toolkit will not scale well beyond around 50k reads/sample. See UMIC-seq pipelines for efficient UMI-gene dictionary generation.

### Profile inserts

- Prepare `data/profile_inserts/sample_probes.csv` with `upstream` and `downstream` columns.
- Run:
  ```bash
  uht-tooling profile-inserts \
    --probes-csv data/profile_inserts/sample_probes.csv \
    --fastq data/profile_inserts/*.fastq.gz \
    --output-dir results/profile_inserts/
  ```
- Outputs: extracted insert FASTA files, QC plots, metrics, and logs. Adjust fuzzy matching strictness via `--min-ratio`.

### EP library profiler (no UMIs)

- Inputs:
  - `data/ep-library-profile/region_of_interest.fasta`
  - `data/ep-library-profile/plasmid.fasta`
  - FASTQ inputs (`--fastq` accepts multiple files)
- Run:
  ```bash
  uht-tooling ep-library-profile \
    --region-fasta data/ep-library-profile/region_of_interest.fasta \
    --plasmid-fasta data/ep-library-profile/plasmid.fasta \
    --fastq data/ep-library-profile/*.fastq.gz \
    --output-dir results/ep-library-profile/
  ```

**Output structure**

Each sample produces an organized output directory:

```
sample_name/
├── KEY_FINDINGS.txt              # Lay-user executive summary
├── summary_panels.png/pdf        # Main visualization
├── aa_mutation_consensus.txt     # Consensus estimate details
├── run.log                       # Analysis log
└── detailed/                     # Technical outputs
    ├── methodology_notes.txt     # Documents which lambda drives what
    ├── lambda_comparison.csv     # Side-by-side lambda comparison
    ├── gene_mismatch_rates.csv
    ├── base_distribution.csv
    ├── aa_substitutions.csv
    ├── plasmid_coverage.csv
    ├── aa_mutation_distribution.csv
    ├── comprehensive_qc_data.csv
    ├── simple_qc_data.csv
    └── qc_plots/                 # QC visualizations
        ├── qc_plot_*.png
        ├── comprehensive_qc_analysis.png
        ├── error_analysis.png
        └── qc_mutation_rate_vs_quality.png/csv
```

**Lambda estimates: which to use**

The profiler calculates lambda (mutations per gene copy) via two methods:

| Method | Formula | Error Quantified? | Used For |
|--------|---------|-------------------|----------|
| Simple | `(hit_rate - bg_rate) × seq_len` | No | KDE plot, Monte Carlo simulation |
| Consensus | Precision-weighted average across Q-scores | Yes | Recommended for reporting |

- **For publication/reporting**: Use the consensus value from `KEY_FINDINGS.txt` or `aa_mutation_consensus.txt`.
- **For understanding distribution shape**: See the KDE plot in `summary_panels.png` (note: uses simple lambda).
- **For detailed error analysis**: See `detailed/comprehensive_qc_data.csv`.

The `KEY_FINDINGS.txt` file provides a plain-language summary including:
- Expected AA mutations per gene copy
- Poisson-based interpretation (% wild-type, % 1 mutation, % 2+ mutations)
- Quality assessment (GOOD/ACCEPTABLE/LOW COVERAGE)

**How the mutation rate and AA expectations are derived**

1. Reads are aligned to both the region of interest and the full plasmid. Mismatches in the region define the "target" rate; mismatches elsewhere provide the background.
2. The per-base background rate is subtracted from the target rate to yield a net nucleotide mutation rate, and the standard deviation reflects binomial sampling and quality-score uncertainty.
3. The net rate is multiplied by the CDS length to estimate λ_bp (mutations per copy). Monte Carlo simulations then flip random bases, translate the mutated CDS, and count amino-acid differences across 1,000 trials—these drive the AA mutation mean/variance that appear in the panel plot.
4. If multiple Q-score thresholds are analysed, the CLI aggregates them via a precision-weighted consensus (1 / standard deviation weighting) after filtering out thresholds with insufficient coverage; the consensus value is written to `aa_mutation_consensus.txt` and plotted as a horizontal guide.

---

## GUI quick start (optional)

The Gradio GUI wraps the same workflows with upload widgets and result previews. Launch it directly:

```bash
python -m uht_tooling.workflows.gui
```

Key points:
- The server binds to `http://127.0.0.1:7860` by default and falls back to an available port if 7860 is busy. Copy http://127.0.0.1:7860 into your browser to interface with the GUI.
- Temporary working directories are created under the system temp folder and cleaned automatically.
- Output archives (ZIP files) mirror the directory structure produced by the CLI.

### Tabs and capabilities

1. **Nextera XT** – forward/reverse primer inputs with CSV preview.
2. **SLIM** – template/context FASTA text areas plus mutation list (supports library codons like `R57:NNK`).
3. **KLD** – inverse-PCR primer design using the same mutation list format (including library codons like `R57:NNK`).
4. **Gibson** – multi-mutation support using `+` syntax.
5. **Mutation Caller** – upload FASTQ and template FASTA, then enter flanks and gene length bounds inline.
6. **UMI Hunter** – long-read UMI clustering with flank entry, UMI length bounds, mutation threshold, and minimum cluster size.
7. **Profile Inserts** – interactive probe table plus multiple FASTQ uploads with adjustable fuzzy-match ratio.
8. **EP Library Profile** – FASTQ uploads plus plasmid and region FASTA inputs.

### Workflow tips

- For large FASTQ datasets, the CLI remains the most efficient option (especially for automation or batch processing).
- Use the command-line flag `--share` in `python -m uht_tooling.workflows.gui` if you need to expose the GUI outside localhost.

### Troubleshooting

- **Port already bound:** the launcher automatically selects the next free port and logs the chosen URL.
- **Missing dependency:** ensure you installed with `pip install "uht-tooling[gui]"`.
- **Stopping the server:** press `Ctrl+C` in the terminal session running the GUI.

---

## Logging

Every workflow configures logging to the destination output directory. Inspect `run.log` for command echoes, parameter choices, and any warnings produced during execution. When providing bug reports, include this log file along with input metadata to streamline triage.

---

## Roadmap

- Replace deprecated Biopython command-line wrappers with native subprocess implementations.
- Expand CLI coverage to any remaining legacy scripts that are still invoked via `make`.
- Add documentation for automation pipelines and integrate continuous integration tests.

Contributions in the form of bug reports, pull requests, or feature suggestions are welcome. File issues on GitHub with clear reproduction steps and sample data when possible.

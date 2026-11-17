# bd-duncani-annotation

Custom scripts for the transcript-defined annotation of the *Babesia duncani* WA1 genome (v3.10), reported in the Data Descriptor *"A transcript-defined annotation of the Babesia duncani genome with measured UTRs and a structure layer"*.

This repository contains **code only**. Assemblies and annotations are deposited at GenBank; all data tables, structure models and derived files are at [Zenodo (DOI 10.5281/zenodo.24101801)](https://zenodo.org/records/24101801).

## What this repository covers

| Directory | Contents |
|---|---|
| `01_assembly/` | Genome comparison against published assemblies (synteny input). |
| `02_annotation/` | GTF construction and curation: UTR reconstruction, TSS/PAS, lncRNA/pseudogene classification, ID renumbering, product-name sanitisation, full-length read validation. |
| `03_structure/` | Boltz-2 structure prediction, Foldseek searches, PDB/sequence consistency checks, pLDDT statistics. |
| `04_validation/` | BUSCO, protein-level QC, and audit scripts that verify manuscript numbers and table gene IDs against the final GTF. |
| `05_figures/` | Scripts that generate Figures 1–5 of the Data Descriptor. |
| `docs/` | Pipeline overview and third-party software versions. |

Scripts for the differential-expression, gene-family, docking and Tpr analyses are **not** included here — they belong to a separate study.

## Dependencies

Most scripts are Python 3 and use `numpy`, `pandas`, `scipy`, `matplotlib`, `biopython`; a few use `Biopython` only. Create the environment with:

```bash
conda env create -f environment.yml
conda activate bd-annotation
```

Third-party command-line tools used by the pipeline (BUSCO, minimap2, BLAST+, hmmscan, Boltz-2, Foldseek, MAFFT, …) are listed with exact versions in [`docs/versions.md`](docs/versions.md). They are invoked by the relevant scripts but are not Python dependencies.

## Reproducing

Scripts are organised by pipeline stage. Each script is self-contained and takes its inputs via command-line arguments or environment variables. Run any script with `--help` to see its options.

Path variables used in the scripts and their meaning:

| Variable / placeholder | Meaning |
|---|---|
| `GEN_ROOT` | Root of the local project tree (contains `V2/`, `annotation_qc/`, `papers/`). |
| `WORK_DIR` | Scratch/working directory for intermediate results. |
| `PAPER_DIR` | Directory holding the manuscript tables used for figures. |
| `TABLE_DIR` | Directory for output tables. |
| `NCRNA_DIR` | Directory holding tRNAscan-SE / barrnap ncRNA predictions. |

A typical figure regeneration looks like:

```bash
export GEN_ROOT=/path/to/your/project
cd 05_figures
python fig_paper1_resource.py    --root "$GEN_ROOT" --out ./figs
python fig_paper1_structure.py   --root "$GEN_ROOT" --out ./figs
python fig_utr_promoter.py       --paper-dir "$PAPER_DIR"   # writes $PAPER_DIR/figs
python fig_lncrna.py             --paper-dir "$PAPER_DIR"   # writes $PAPER_DIR/figs
```

Note that `fig_utr_promoter.py` and `fig_lncrna.py` also expect the submitted
tables under `$PAPER_DIR/tables` and the validation tables under
`$PAPER_DIR/../Submission_Package/02_data/tables`.

The scripts that produce the final annotation are staged in `02_annotation/`; [`docs/pipeline.md`](docs/pipeline.md) describes the order and the input/output of each step. Because the pipeline starts from raw Nanopore/Illumina reads, a full from-scratch reproduction requires the read data from SRA; the committed tables and the final GTF at Zenodo let you re-run any individual validation or figure step directly.

## Citation

If you use this code, please cite the Data Descriptor (manuscript DOI to be added on acceptance) and the archived release:

> Wang S. *et al.* A transcript-defined annotation of the *Babesia duncani* genome with measured UTRs and a structure layer. (manuscript DOI pending).

> Zenodo archive: <https://zenodo.org/records/24101801> (DOI 10.5281/zenodo.24101801), release v3.10.

## Licence

MIT — see [LICENSE](LICENSE).

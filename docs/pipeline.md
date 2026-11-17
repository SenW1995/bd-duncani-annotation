# Pipeline overview

This Data Descriptor reports a transcript-defined annotation of the *Babesia duncani* WA1 genome. The pipeline starts from raw long/short reads and ends at the v3.10 annotation, its validation tables and Figures 1–5. The stages below map to the directories in this repository.

## 1. Assembly curation (`01_assembly/`)

- Input: hifiasm contigs, Nanopore/Illumina reads, Hi-C reads.
- Third-party steps (not custom code): hifiasm assembly → nextPolish2 polishing → Juicer/3D-DNA or YaHS Hi-C scaffolding → removal of apicoplast/mitochondrial contigs → Merqury QV.
- `srv_17_genome_compare.sh` prepares the pairwise alignments used for the synteny panel of Fig. 1 against the Singh 2023 and Wang 2022 assemblies.

## 2. Annotation construction (`02_annotation/`)

Run in roughly this order:

1. `t3_tss.py` — transcription start sites and unannotated transcribed loci.
2. `srv_05_utr.py` — rebuild 5′/3′ UTRs from full-length Nanopore cDNA; emits `five_prime_utr` / `three_prime_utr` GTF features.
3. `srv_11_polya.py` — polyadenylation sites (PAS) and alternative polyadenylation.
4. `srv_20_lncrna_final.py` — final lncRNA set and pseudogene classification.
5. `srv_55_merge_ncrna.py` — merge tRNAscan-SE and barrnap ncRNA calls into the annotation → **v3.10 GTF**.
6. `srv_29_renumber_ids.py` — unified `BdWA1_###g#####` gene numbering.
7. `srv_73_fix_frameshift.py` — convert coding models with unresolved frameshifts to pseudogenes.
8. `srv_74_sanitize_products.py` — normalise `product` qualifiers for GenBank submission.
9. `fl_validate.py` — full-length read support; exports the final intron/support tables.
10. `validate_gtf.py`, `translate_cds_qc.py` — CDS re-extraction, translation and integrity checks.
11. `counts_to_tpm.py`, `summarize_resource_numbers.py` — expression normalisation and headline resource counts.

## 3. Structure layer (`03_structure/`)

- `srv_07b_boltz.py`, `srv_22_fold_all.py` — Boltz-2 structure prediction (whole proteome, one representative per multi-copy family).
- `srv_35_foldseek_api.py` — Foldseek searches for unknown-function proteins.
- `srv_43_verify_pdb_seq.py`, `srv_44_verify_v38_protein.py` — verify each PDB matches its annotated protein sequence.
- `srv_45_apply_struct_annot.py` — transfer trusted structural neighbours into product annotations.
- `srv_56_all_plddt.py` — per-model mean pLDDT over the full model set.

## 4. Validation (`04_validation/`)

- `busco_compare.sh` — BUSCO completeness (genome and protein modes, apicomplexa_odb10).
- `protein_level_qc.py` — protein-level QC against 11 apicomplexan reference proteomes plus human.
- `audit_paper_consistency.py` — recomputes every number quoted in the manuscript from the final GTF.
- `audit_table_ids.py` — checks that every gene ID in the submitted tables exists in the final GTF.

## 5. Figures (`05_figures/`)

| Script | Figure |
|---|---|
| `fig_paper1_resource.py` | Fig. 1 (AT content, synteny) and Fig. 2 (gene structure / transcript support) |
| `fig_utr_promoter.py` | Fig. 3 (UTRs, TSS promoters, PAS) |
| `fig_lncrna.py` | Fig. 4 (lncRNA layer and cis relationships) |
| `fig_paper1_structure.py` | Fig. 5 (structure layer and structural annotations) |
| `pub_fig_synteny.py` | three-track synteny panel used in Fig. 1 |
| `pub_style.py` | shared plotting style and panel/alignment helpers |

`pub_style.py` optionally uses an external panel-alignment helper for gate checks; set `SKILL_SCRIPTS` to its directory to enable. If unset, the gate prints a warning and figure generation proceeds.

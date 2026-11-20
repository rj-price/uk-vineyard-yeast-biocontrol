# Analysis code: UK vineyard yeasts as biocontrol agents

Analysis and figure-generation code for the study of UK vineyard yeast isolates
screened for biological control activity against *Botrytis cinerea*,
*Fusarium graminearum* and *Fusarium oxysporum*, combining phenotypic screening
with whole-genome sequencing, assembly, species identification and comparative
genomics.

This repository is deposited to satisfy the *Access Microbiology* open-data policy.
It contains the code used to produce the results, tables and figures in the
manuscript. Software packages and versions are given in the manuscript Methods; the
main ones are also summarised below.

- **Manuscript:** Darko & Price, *Genome mining of UK vineyard yeasts partly assigns
  biocontrol mechanism between genera but not strain-level antagonism*
  (submitted to *Access Microbiology*).
- **Sequence data:** NCBI BioProject **PRJNA1502812** (55 genome assemblies and raw
  Illumina reads). BioSample accessions are listed in Table S1.

## What is and isn't here

This is a **code** repository. It contains the analysis scripts plus the small input
manifests they need (sample sheets, reference-genome accession lists, the assembly
manifest). It does **not** contain:

- raw or assembled sequence data, or reads (these are at NCBI, PRJNA1502812);
- downloaded public reference genomes (retrieve them from the accession lists in
  `phylo_refs/genomes.tsv`, `species_id/ref_accessions.tsv`,
  `species_id/metsch_clade/acc.txt` and `global_context/` with NCBI `datasets`);
- large intermediate outputs, tool result files, SLURM logs and compiled binaries;
- the phenotype spreadsheets and processed data tables (these accompany the paper
  as supplementary files).

The scripts are a record of the analysis. Paths inside them point at the original
working tree on a SLURM HPC cluster, so they will not run unedited from a fresh
clone.

## Repository layout

Shell scripts (`*.sh`) are SLURM batch jobs; Python scripts (`*.py`) do the
downstream analysis, statistics, figures and tables.

| Directory | Stage |
|-----------|-------|
| `reassembly/` | Platanus-allee assembly of the heterozygous *H. uvarum* and *P. kluyveri* genomes, reads downsampled to ~70x (`05_platanus.sh`). |
| `phylo_refs/` | BUSCO on isolates and reference genomes, then the BUSCO supermatrix and IQ-TREE phylogeny (Fig 2A). |
| `species_id/` | Genome-based species identification: fastANI against type and reference genomes (`01`, and the expanded panel in `04` used for Fig 2B and Fig S1), the *M. pulcherrima* clade comparison (`02`), and LSU D1/D2 barcodes for taxa without a close reference genome (`03`). |
| `mining/` | funannotate annotation, dbCAN, MEROPS, antiSMASH (`01`, `02`) and the SignalP secretome (`03`). `killer_subjects.tsv` is the isolate-to-assembly manifest. |
| `interpro/` | InterProScan annotation of the mined genomes (GO terms for the clonal-variant tests). |
| `allelic/` | Within-species analyses: within-species fastANI (`02`), read-based variant calling in the two near-clonal *Rhodotorula* groups (`06`, `07`; BWA and bcftools), the callable genome (`10`), and accessory genes (`08`; CD-HIT and miniprot). |
| `global_context/` | Kent isolates against global *H. uvarum* and *Rhodotorula* genomes (fastANI and the *H. uvarum* BUSCO tree; Fig 5). Repeat deposits of the same strain are listed once, per `meta/duplicate_assemblies.tsv`. |
| `scripts/` | Python pipeline: master table, statistics, species calls, gene inventories, clonal variants, all figures and tables. |
| `exploratory/` | Analyses that were run during the project but are not reported in the paper (see below). |

## Scripts behind each figure and table

| Output | Script |
|--------|--------|
| Fig 1 | `scripts/19_fig1_screen.py` |
| Fig 2 | `scripts/20_fig2_phylo.py` |
| Fig 3 | `scripts/21_fig3_mechanisms.py`; iron statistics in `scripts/39_iron_stats.py` |
| Fig 4 | `scripts/22_fig4_genomic.py` |
| Fig 5 | `scripts/25_fig5_global.py` |
| Fig 6 | `scripts/24_fig6_prioritisation.py`; composite score in `scripts/02_stats_figures.py`, robustness in `scripts/37_score_sensitivity.py` |
| Fig S1 | `scripts/28b_figS1_sporo_ani.py` |
| Fig S2 | `scripts/23_figS2_clonal_ani.py` |
| Fig S3 | `scripts/26_figS3_biocontrol_global.py` |
| Tables S1 to S9 | `scripts/27_supp_tables.py` |
| Clonal variants (121 sites, Table S8) | `scripts/38_clonal_reads.py`, `16_clonal_variants.py`, `18_interpro_enrichment.py` |
| Accessory genes | `scripts/36_accessory.py` |

`scripts/16_clonal_variants.py` and `18_interpro_enrichment.py` default to the UY1
group; `CLONAL_REF=TY28 CLONAL_QRY=UY2+UY3` runs the second group.

Helper modules: `scripts/_palette.py` (shared species colours), `scripts/_figstyle.py`
(matplotlib style), `scripts/_mk_*_sheet.py` (sample-sheet builders).
`scripts/30_ncbi_upload.py` to `34_fill_genome_template.py` prepared the NCBI
deposition (assembly renaming and cleaning, metadata); they are not analyses.

## Software and environments

Tools were run inside per-task conda environments (`conda run -n <env> <command>`).

| Environment | Tools |
|-------------|-------|
| `data_science` | Python 3.13, pandas, NumPy, SciPy, statsmodels, Matplotlib, Biopython |
| `funannotate` | funannotate 1.8.17 (gene prediction and annotation), openpyxl |
| `busco` / `BUSCO_phylogenomics` | BUSCO 5 (`fungi_odb10`), supermatrix building; IQ-TREE 2 |
| `fastani` | fastANI 1.34 |
| `extract_its` | barrnap (rRNA gene extraction for D1/D2 barcodes) |
| `blast` | blastn (D1/D2 barcodes against NCBI nt) |
| `variant_calling` | BWA, SAMtools, bcftools (clonal variants) |
| `cdhit`, `evm` | CD-HIT, miniprot (accessory genes) |
| `run_dbcan` | dbCAN 5.0 |
| `antismash` | antiSMASH 8.0.1 |
| `interproscan` | InterProScan 5.74-105.0 (`-goterms -iprlookup`) |

Other tools called by the SLURM scripts: Platanus-allee 2.4.0, fastp, seqtk and
SignalP 5.0b.

## Run order (Python pipeline)

The upstream SLURM stages (assembly, annotation, mining, phylogeny, ANI, variant
calling) run first; `scripts/` consume their outputs.

```bash
conda run -n funannotate  python scripts/00_dump_xlsx.py        # phenotype sheets -> CSV
conda run -n data_science python scripts/01_build_master.py     # isolate master table
conda run -n data_science python scripts/02_stats_figures.py    # phenotype summaries, composite score
conda run -n data_science python scripts/07_platanus_comparison.py
conda run -n data_science python scripts/08_mining_summary.py   # gene inventory
conda run -n data_science python scripts/10_secretome.py
conda run -n data_science python scripts/14_species_final.py    # species calls
conda run -n data_science python scripts/38_clonal_reads.py     # read-validated clonal variants
conda run -n data_science python scripts/16_clonal_variants.py
conda run -n data_science python scripts/18_interpro_enrichment.py
conda run -n data_science python scripts/36_accessory.py
conda run -n data_science python scripts/37_score_sensitivity.py
conda run -n data_science python scripts/39_iron_stats.py
# figures: 19-26, 28b; tables: 27
```

## Exploratory analyses

`exploratory/` keeps work that was run but did not make it into the paper, so the
record is complete. None of it is needed to reproduce the reported results.

| Path | What it was |
|------|-------------|
| `exploratory/mining/04_killer_tblastn.sh`, `killer_toxins.faa`, `exploratory/scripts/12_killer_summary.py` | Killer-toxin tblastn search |
| `exploratory/allelic/01_orthofinder.sh` | OrthoFinder orthogroups |
| `exploratory/allelic/03_clonal_pair.sh`, `05_clonal_ref.sh` | Assembly-based dnadiff comparison of UY1 and TY27, replaced by the read-based calls. `scripts/38_clonal_reads.py` still intersects its output with the read calls as a check |
| `exploratory/allelic/04_allelic.sh`, `exploratory/scripts/15_allelic.py` | Allelic and upstream variants in *H. uvarum* |
| `exploratory/allelic/09_rhodo8_ani.sh` | Eight-isolate *Rhodotorula* fastANI; Fig S2 uses the all-versus-all run in `global_context/` |
| `exploratory/mapping_qc/`, `exploratory/scripts/05_mapping_qc_summary.py` | Read mapping to diagnose the fragmented ascomycete SPAdes assemblies |
| `exploratory/reassembly/`, `exploratory/scripts/06_reassembly_summary.py`, `28_figS3_assembly.py` | SPAdes re-assembly, GenomeScope k-mer profiles and an assembly-rescue figure |
| `exploratory/scripts/03_phylogeny.py`, `04_genomics_qc.py`, `09_genotype_phenotype.py`, `11_comparative.py`, `13_species_id.py` | Early versions of the tree, assembly QC (incl. kraken2), genotype-phenotype tables and species calls |
| `exploratory/scripts/17_regulatory.py` | PFAM-based regulatory test, replaced by the InterPro/GO version (`18`) |
| `exploratory/scripts/35_phylo_signal.py` | Parsimony test of phylogenetic signal in *H. uvarum* |
| `exploratory/scripts/29_build_docx.py` | Draft document builder |
| `exploratory/global_context/probe_*.py` | Ad hoc inspection of ANI tables, trees and metadata |

## Isolate naming key

Three code systems refer to the same isolates:

| System | Example | Where used |
|--------|---------|------------|
| Assay code (paper) | `UY1`, `TY26`, `RY3` | phenotype data, manuscript |
| Genome/culture code | `CUY1`, `CTY26`, `RY3` | culture collection |
| Sequencing ID | `Y1_1`, `Y2_9`, `Y3_5` | raw reads, tree tips |

The genome code is the assay code with a leading culture `C` (`CUY` to `UY`, `CTY` to
`TY`, `RY` unchanged); the sequencing ID maps via the master table
(`scripts/01_build_master.py`) and Table S5.

## Citation

If you use this code, please cite the paper and the archived release (see
[`CITATION.cff`](CITATION.cff)).

## Licence

Released under the MIT Licence; see [`LICENSE`](LICENSE).

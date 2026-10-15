# Analysis code — UK vineyard yeasts as biocontrol agents

Analysis and figure-generation code for the study of UK vineyard yeast isolates
screened for biological control activity against *Botrytis cinerea*,
*Fusarium graminearum* and *Fusarium oxysporum*, combining phenotypic screening
with whole-genome sequencing, assembly, species identification and comparative
genomics.

This repository is deposited to satisfy the *Access Microbiology* open-data policy:
it contains **all of the code** used to produce the results, tables and figures in
the manuscript. Software packages and versions are described in the manuscript
Methods; the key ones are also summarised below.

- **Manuscript:** Darko & Price, *Genome mining of UK vineyard yeasts partly assigns biocontrol mechanism between genera but does not predict antagonist strength within a species*
  (submitted to *Access Microbiology*).
- **Sequence data:** NCBI BioProject **PRJNA1502812** (55 genome assemblies + raw
  Illumina reads). Genome accessions will be added on release.

## What is and isn't here

This is a **code** repository. It contains every analysis script plus the small
input manifests they need (sample sheets, reference-genome accession lists, the
curated killer-toxin protein panel). It deliberately does **not** contain:

- raw or assembled sequence data, or reads — these are at NCBI (PRJNA1502812);
- downloaded public reference genomes — retrieve them from the accession lists in
  `phylo_refs/genomes.tsv`, `species_id/ref_accessions.tsv` and
  `global_context/` using NCBI `datasets`;
- large intermediate outputs, tool result files, SLURM logs and compiled binaries;
- the phenotype spreadsheets and processed data tables (these accompany the paper
  as supplementary files).

The scripts are provided as a faithful record of the analysis. Paths inside them
point at the original working tree on a SLURM HPC cluster, so they are not intended
to run unedited from a fresh clone — they document exactly how each result was
produced.

## Repository layout

The layout mirrors the analysis pipeline. Shell scripts (`*.sh`) are SLURM batch
jobs run on the cluster; Python scripts (`*.py`) do the downstream analysis, stats
and figures.

| Directory | Stage |
|-----------|-------|
| `scripts/` | Main Python pipeline: build the isolate master table, phenotype stats, QC, species ID, comparative genomics, and all manuscript figures/tables (numbered run order). |
| `phylo_refs/` | Reference-genome BUSCO + supermatrix phylogeny (`genomes.tsv` = reference set). |
| `reassembly/` | Ascomycete re-assembly: k-mer/GenomeScope2 heterozygosity analysis and Platanus-allee haplotype-aware assembly. |
| `mapping_qc/` | fastp QC + read-mapping back to assemblies to diagnose the ascomycete assembly failure. |
| `species_id/` | Genome-based species identification (fastANI vs type genomes + LSU D1/D2 barcode BLAST), and TY26/TY30 against the *Metschnikowia pulcherrima* clade (`02_metsch_clade.sh`, accessions in `metsch_clade/acc.txt`). |
| `global_context/` | Placing the UK isolates against global *Hanseniaspora uvarum* / *Rhodotorula* references (ANI + tree). Re-deposits of the same strain are listed once, per `meta/duplicate_assemblies.tsv`. |
| `mining/` | Genome annotation (funannotate) and biocontrol-gene mining: CAZymes (dbCAN), secondary-metabolite clusters (antiSMASH), SignalP secretome, killer-toxin tblastn (`killer_toxins.faa` = curated panel). |
| `interpro/` | InterProScan functional re-annotation of the mined genomes. |
| `allelic/` | Within-species comparative genomics: OrthoFinder, within-clade ANI, eight-isolate *Rhodotorula* ANI (`09`), read-based variant calling within the two near-clonal *Rhodotorula* groups (`06`, `07`; BWA + bcftools) and accessory genes (`08`; CD-HIT + miniprot). The assembly-based dnadiff comparison (`03`, `05`) was superseded by the read-based calls. |

Helper modules: `scripts/_palette.py` (shared species→colour map for every figure),
`scripts/_figstyle.py` (matplotlib style), `scripts/_mk_*_sheet.py` (spreadsheet
builders). `scripts/30_ncbi_upload.py` and `31`–`34` are NCBI data-deposition helpers (assembly
cleaning and metadata generation), not analyses.

## Software and environments

Tools were run inside per-task conda environments (`conda run -n <env> <command>`)
on the cluster. Principal software (full version list in the manuscript Methods):

| Environment | Tools |
|-------------|-------|
| `data_science` | Python 3, pandas, numpy, scipy, statsmodels, matplotlib, biopython |
| `funannotate` | funannotate 1.8.17 (gene prediction + annotation), openpyxl |
| `busco` / `BUSCO_phylogenomics` | BUSCO (`fungi_odb10`), supermatrix building; ML tree with IQ-TREE 2 |
| `fastani` | fastANI (genome-based species ID and within-clade ANI) |
| `variant_calling` | BWA, SAMtools, bcftools (read-based clonal variants) |
| `cdhit`, `evm` | CD-HIT, miniprot (accessory genes) |
| `mummer` | MUMmer `dnadiff` (initial clonal-pair comparison, superseded) |
| `orthofinder` | OrthoFinder (orthogroup inference) |
| `run_dbcan` | dbCAN 5.0 (CAZyme annotation) |
| `antismash` | antiSMASH 8.0.1 (secondary-metabolite clusters) |
| `interproscan` | InterProScan 5.74-105.0 (`-goterms -iprlookup`) |
| `blast` | tblastn (killer-toxin search), blastn (D1/D2 barcodes) |

Other tools invoked by the SLURM scripts: SPAdes 3.15.5 and Platanus-allee 2.4.0
(assembly), fastp, BWA and samtools (mapping QC), Jellyfish + GenomeScope2 (k-mer
heterozygosity), barrnap (rRNA extraction), SignalP 5.0b (secretome).

## Run order (Python pipeline)

The `scripts/` are numbered in dependency order. Representative sequence:

```bash
conda run -n funannotate  python scripts/00_dump_xlsx.py       # phenotype sheets -> CSV
conda run -n data_science python scripts/01_build_master.py    # isolate master table
conda run -n data_science python scripts/02_stats_figures.py   # phenotype stats
# ... 03–18: QC, phylogeny, species ID, mining, secretome, comparative genomics ...
# 38_clonal_reads.py (read-validated clonal variants) runs before 16–18;
# 35–37: lineage signal, accessory genes, prioritisation-score sensitivity
conda run -n data_science python scripts/19_fig1_screen.py     # manuscript figures
# ... 20–28: remaining main + supplementary figures ...
conda run -n data_science python scripts/27_supp_tables.py     # supplementary tables
conda run -n funannotate  python scripts/29_build_docx.py      # assemble document
```

The upstream SLURM stages (assembly, annotation, mining, phylogeny, ANI) must be run
first; the `scripts/` consume their outputs.

## Isolate naming key

Three code systems refer to the same physical isolates:

| System | Example | Where used |
|--------|---------|------------|
| Assay code (paper) | `UY1`, `TY26`, `RY3` | phenotype data, manuscript |
| Genome/culture code | `CUY1`, `CTY26`, `RY3` | culture collection |
| Sequencing ID | `Y1_1`, `Y2_9`, `Y3_5` | raw FASTA / tree tips |

The genome code is the assay code with a leading culture-`C` (`CUY→UY`, `CTY→TY`,
`RY→RY`); the sequencing ID maps via the master table (`scripts/01_build_master.py`).

## Licence

Released under the MIT Licence — see [`LICENSE`](LICENSE).

# Treg chromatin accessibility in the tumour microenvironment

Analysis code for:

> **Regulatory T cells in the tumor microenvironment display a unique chromatin
> accessibility profile**
> Dadey R.E., Cohen J.A., Cui J., *et al.* *ImmunoHorizons* (2025).
> doi:[10.1093/immhor/vlae014](https://doi.org/10.1093/immhor/vlae014)

This repository contains the ATAC-seq and RNA-seq analysis code behind
**Figures 1–4** and the sequencing quality-control supplements. Wet-lab
procedures (mouse work, sorting, library preparation) are described in the
paper and have no code here.

---

## What the analysis does

Tregs and CD4<sup>+</sup>Foxp3<sup>−</sup> T cells were sorted from B16 and MC38
tumours and from non-draining lymph nodes, profiled by low-input ATAC-seq and
RNA-seq, and compared along three axes:

1. **Location vs time** — is intratumoural Treg chromatin accessibility driven
   more by tissue or by tumour timepoint? (Fig 1)
2. **Tumour vs lymph node** — which regions open and close, and do they carry a
   Treg effector signature? (Fig 2)
3. **Tumour-specific signature** — subtracting skin, non-lymphoid and
   resting/activated Treg profiles to isolate peaks unique to tumour Tregs.
   (Fig 3)
4. **Effect of immunotherapy** — does anti-PD-1 treatment alter the Treg
   chromatin or transcriptional profile? (Fig 4)

---

## Repository layout

```
1.Data/     input data (not distributed — see "Data availability")
2.Script/   analysis scripts, numbered in figure order
3.Output/   tables and figures produced by the scripts (not distributed)
```

Scripts expect the **repository root** as the R working directory:

```r
setwd("/path/to/treg-tme-chromatin")
source("2.Script/01_samplesheet_generation.R")
```

Each script carries a header giving its figure panel, inputs, outputs, the
upstream script it depends on, and its original filename in the lab archive.

---

## Figure → script map

| Figure | Analysis | Script |
| --- | --- | --- |
| — | nf-core sample sheets for the ATAC-seq and RNA-seq runs | [`01_samplesheet_generation.R`](2.Script/01_samplesheet_generation.R) |
| — | Merge sample metadata across ATAC-seq batches | [`02_metadata_integration.R`](2.Script/02_metadata_integration.R) |
| **1B–1E** | Differential accessibility and PCA, tumour vs LN; replicate correlation | [`03_atac_deseq_tumor_vs_ln.R`](2.Script/03_atac_deseq_tumor_vs_ln.R) |
| **2B** | ChIPseeker peak annotation and genomic feature distribution | [`04_atac_peak_annotation.R`](2.Script/04_atac_peak_annotation.R) |
| **S1A–B** | TSS enrichment scores and ATAC-seq mapping statistics | [`05_atac_tss_enrichment_qc.R`](2.Script/05_atac_tss_enrichment_qc.R) |
| Suppl table | STAR mapping statistics for RNA-seq libraries | [`06_rnaseq_mapping_stats.R`](2.Script/06_rnaseq_mapping_stats.R) |
| **2A, 2C–2J** | Opening/closing peak counts, GSEA vs Treg signatures, KEGG enrichment | [`07_enrichment_tumor_vs_ln.R`](2.Script/07_enrichment_tumor_vs_ln.R) |
| **2** | Accessibility of Treg signature genes around the TSS | [`08_treg_signature_tss_viz.R`](2.Script/08_treg_signature_tss_viz.R) |
| **3A** | PCA with public non-lymphoid and resting/activated Tregs | [`09_atac_deseq_public_datasets.R`](2.Script/09_atac_deseq_public_datasets.R) |
| **3B** | Differential accessibility and PCA of skin Tregs | [`10_atac_deseq_skin_treg.R`](2.Script/10_atac_deseq_skin_treg.R) |
| **3B** | Joint PCA of skin, B16 and MC38 Tregs | [`11_atac_deseq_skin_b16_mc38.R`](2.Script/11_atac_deseq_skin_b16_mc38.R) |
| **3C, 3H** | Subtraction method isolating tumour-specific peaks | [`12_subtraction_tumor_specific.R`](2.Script/12_subtraction_tumor_specific.R) |
| **3H** | Venn comparison of skin vs tumour-specific peaks and genes | [`13_compare_skin_vs_tumor.R`](2.Script/13_compare_skin_vs_tumor.R) |
| **3** | KEGG over-representation for skin-only, tumour-only and shared genes | [`14_enrichment_skin_tumor_treg.R`](2.Script/14_enrichment_skin_tumor_treg.R) |
| **3D–3G** | GSEA of the tumour-specific list against the tissue Treg signature | [`15_gsea_subtraction_list.R`](2.Script/15_gsea_subtraction_list.R) |
| **3J** | KEGG and GO over-representation of the tumour-specific gene list | [`16_enrichment_tumor_vs_skin.R`](2.Script/16_enrichment_tumor_vs_skin.R) |
| **3I–3J** | GO bar plot and fold-change visualisation | [`17_go_barplot_and_deg_viz.R`](2.Script/17_go_barplot_and_deg_viz.R) |
| **4B, 4D, 4F** | ATAC-seq after anti-PD-1 or isotype treatment | [`18_atac_deseq_antipd1.R`](2.Script/18_atac_deseq_antipd1.R) |
| **4B, 4D, 4F** | Same, after ComBat-seq batch correction | [`19_atac_deseq_antipd1_batchcor.R`](2.Script/19_atac_deseq_antipd1_batchcor.R) |
| **4C, 4E, 4G** | RNA-seq after anti-PD-1 or isotype treatment | [`20_rnaseq_deseq_antipd1.R`](2.Script/20_rnaseq_deseq_antipd1.R) |
| **4C, 4E, 4G** | Same, after ComBat-seq batch correction | [`21_rnaseq_deseq_antipd1_batchcor.R`](2.Script/21_rnaseq_deseq_antipd1_batchcor.R) |

Run `01` → `02` → `03` → `04` first; `04` writes the annotated peak tables that
every enrichment script reads, and `12` writes the tumour-specific peak list
that Figures 3D–3J depend on.

---

## Data availability

**No sequencing data, count matrices or intermediate objects are included here.**

| Dataset | Source |
| --- | --- |
| ATAC-seq and RNA-seq generated for this study | GEO [**GSE246703**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246703) |
| Non-lymphoid tissue Treg ATAC-seq (public) | GEO [GSE112731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112731) |
| Resting / activated Treg ATAC-seq (public) | GEO [GSE83315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83315) |

### Upstream processing happens before these scripts

The scripts in this repository start from **count matrices**, not FASTQ files.
Read alignment and peak calling were run separately, as described in the paper:

- **ATAC-seq** — BWA-MEM to mm10, peak calling with MACS2 2.1.2
  (`-B -q 0.01 --nomodel -f BAM -g mm`), blacklisted regions removed,
  quantification with `Rsubread::featureCounts` 1.28.1. The scripts expect the
  consensus peak count table
  (`1.Data/ATAC_overall/consensus_peaks.mLb.clN.featureCounts.txt`).
- **RNA-seq** — STAR 2.5.1 to GRCm38, counts via HTSeq-count. The scripts
  expect the per-sample `ReadsPerGene.out.tab` files in `1.Data/RNA_readTabls/`.

The public datasets must be downloaded and put through the same pipeline before
running `09_atac_deseq_public_datasets.R`.

### Reference gene list

Scripts `07`, `08` and `15` read a Treg effector/transcriptional signature table
at `1.Data/Treg_effector_sig.xlsx` (columns: `EntrezID`, `Symbol`, `AveExpr`,
`t`, `P.Value`, `adj.P.Val`). It is a published differential-expression result
used as the GSEA reference set and is not redistributed here.

---

## Software

| Tool | Version |
| --- | --- |
| BWA-MEM | — |
| MACS2 | 2.1.2 |
| Rsubread (`featureCounts`) | 1.28.1 |
| STAR | 2.5.1 |
| HTSeq-count | — |
| DESeq2 | — |
| ChIPseeker | — |
| clusterProfiler | — |
| fgsea / msigdbr | MSigDB KEGG, C2, C5, C7 |
| ataqv | — |
| sva (`ComBat_seq`) | — |

Multiple testing was controlled with Benjamini–Hochberg; FDR ≤ 5 % was used
throughout.

---

## Contact

**Jian Cui** — Senior Bioinformatics Scientist  
Vignali Lab, Department of Immunology  
School of Medicine, University of Pittsburgh  
<cuij@pitt.edu>

Questions about the analysis or this code are welcome.

## License

[MIT](LICENSE)

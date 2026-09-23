# OsMYB80 / TDR anther RNA-seq analysis

R code for the anther RNA-seq analysis of wild-type, *osmyb80* and *tdr* rice (*Oryza sativa* cv. Nipponbare) in

> **[Authors] ([Year]) [Title]. [Journal]. [DOI]**

It reproduces the differential expression analysis, the **HOT (high temperature–OsMYB80–TDR) genes**, and the Fig. 5 heat map. The whole analysis is in a single script, `OsMYB80_TDR_RNAseq_analysis.R`. Its input is the count matrix and sample metadata deposited in GEO (**[GSExxxxxx]**).

## Repository contents

```
OsMYB80_TDR_RNAseq_analysis.R    the analysis (steps 1-5 below)
data/
  HOT_genes_Fig5.csv             the 36 HOT genes shown in Fig. 5 (order used in the heat map)
  high_temperature_repressed_71genes.csv   71 genes of clusters 2 and 15 in Endo et al. (2009)  [to be added]
  GO/goid2gene_BP.txt, godb_BP.txt         rice GO-BP annotation (RAP-DB + Oryzabase), see below
  GO/riceIDtable.csv             RAP-DB ↔ Entrez Gene ID table
input/                           GEO supplementary files (not tracked; downloaded)
results/                         outputs
```

## Input data (GEO [GSExxxxxx])

| file | content |
|---|---|
| `OsMYB80_TDR_anther_Salmon_NumReads.txt.gz` | Salmon-estimated read counts (NumReads); rows = RAP-DB representative transcripts, columns = 56 samples |
| `OsMYB80_TDR_anther_sample_metadata.txt.gz` | genotype, mutant allele, plant ID, auricle distance, panicle length, and column order used in the analysis |

There are 56 libraries: wild type (n = 18), *osmyb80* (n = 10) and *tdr* (n = 28). Each library comes from the anthers of one plant. Place the two files in `input/`, or set `GEO_ACCESSION` so the script downloads them from the GEO FTP site.

### Upstream processing (not part of this repository)

Read processing, mapping, and transcript quantification were performed as described in Kashima et al. (2021):

- read 1 was trimmed with fastp (`--trim_poly_x`, Illumina adapters, `-l 31`);
- reads were mapped with BWA-MEM to the RAP-DB IRGSP-1.0 representative transcript sequences (`IRGSP-1.0_representative_transcript_exon_2022-09-01.fasta`);
- transcripts were quantified with Salmon v1.4.0 in alignment-based mode (`salmon quant -l IU -a <BAM>`).

The `NumReads` column of `quant.sf` is the input of this analysis.

## Requirements

The published results were produced with **R 4.1 / Bioconductor 3.14**:

| package | version |
|---|---|
| Seurat / SeuratObject | 4.2.0 / 4.1.2 |
| TCC / edgeR | 1.34.0 / 3.36.0 |
| clusterProfiler | 4.2.2 |
| AnnotationHub | 3.2.2 |
| ggplot2, patchwork, VennDiagram | – |

**Use these versions.** With newer edgeR (≥ 4.x, which TCC calls internally) the numbers of DEGs change.

A reproducible way to get them is the Bioconductor 3.14 container:

```bash
singularity pull bioc_3.14.sif docker://bioconductor/bioconductor_docker:RELEASE_3_14   # or docker pull
singularity exec --cleanenv --env R_LIBS_USER=$HOME/R/bioc3.14 bioc_3.14.sif R -e '
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/2022-09-29"))
  BiocManager::install(c("TCC", "Seurat", "clusterProfiler", "AnnotationHub", "VennDiagram", "patchwork"),
                       lib = Sys.getenv("R_LIBS_USER"), update = FALSE, ask = FALSE,
                       site_repository = getOption("repos")[["CRAN"]])'
```

## Usage

```bash
# from the repository root
GEO_ACCESSION=GSExxxxxx Rscript OsMYB80_TDR_RNAseq_analysis.R
# or, with the files already in input/
Rscript OsMYB80_TDR_RNAseq_analysis.R
```

The run takes a few minutes. The GO step (AnnotationHub) needs internet access the first time. Set `RUN_GO <- FALSE` in the script to skip it.

## Analysis steps

1. **Input.** Count matrix and metadata are read from the GEO files. Each RAP representative transcript is one feature. Transcripts are not summed per gene: `Os04t0470600-01` → `Os04g0470600`, and further transcripts of the same locus get `.1`, `.2`, …, as named by Seurat.
2. **Seurat.** `NormalizeData` (LogNormalize, scale factor 1 × 10^6), `FindVariableFeatures`, `ScaleData`, `RunPCA`. QC plots are written to `results/QC_seurat_myb80_tdr.pdf`.
3. **Differential expression (TCC).** *osmyb80* vs WT and *tdr* vs WT are analysed separately: `calcNormFactors(iteration = 3)`, then `estimateDE(FDR = 0.05)`. DEGs are the features with `estimatedDEG == 1`, i.e. q < 0.05 at FDR = 0.05. Up or down is decided from the mean TCC-normalised expression in mutant vs WT. `log2FC_mut_vs_WT` is always log2(mutant/WT); TCC's own `m.value` takes the ratio relative to the group of the first sample column, and its raw value is kept as `m.value_TCC_raw`.
4. **HOT genes and Fig. 5.**
   - Genes downregulated in both mutants (**FDR = 0.05, log2FC < −1**) are compared with the 71 high-temperature-repressed genes (Endo et al., 2009). Genes in both sets are the HOT genes.
   - The script writes the overlap counts, the gene lists, the status of each of the 36 HOT genes, and the heat map (`DoHeatmap` of `ScaleData`-scaled expression, all samples).
   - Counts without the fold-change cutoff are also reported.
5. **GO enrichment.**
   - **A:** `enrichGO` with the AnnotationHub OrgDb AH96211 (Entrez IDs), for genes down in both mutants and for genes down in *tdr* only.
   - **B:** HOT genes, using both `enrichGO` and `enricher`. `enricher` uses the RAP-DB (IRGSP-1.0 representative annotation 2021-11-11) + Oryzabase GO-BP annotation, with all quantified loci as background.

## Output (`results/`)

| file | content |
|---|---|
| `DEG_myb80.csv`, `DEG_tdr.csv` | TCC results for all features (normalised means, log2FC, p, q, DEG call, direction) |
| `DEG_summary.csv` | number of up/downregulated features and loci |
| `Fig5_venn_counts.csv`, `gene_sets/` | Venn region counts and gene lists; `Fig5_venn.pdf` when the 71-gene list is present |
| `HOT_genes_status.csv` | TCC statistics of the 36 HOT genes |
| `Fig5_heatmap.pdf` | Fig. 5 heat map |
| `Fig5_heatmap_auricle_subset_ALTERNATIVE.pdf` | variant defined in the original code (subset by auricle distance); not used in the paper |
| `GO/` | GO enrichment tables and dot plots |
| `sessionInfo.txt` | R session information |

Main numbers (features = RAP transcripts):

| comparison | up (FDR = 0.05) | down (FDR = 0.05) | down (FDR = 0.05, log2FC < −1) |
|---|---|---|---|
| *osmyb80* vs WT | 2,106 | 5,291 | 3,274 |
| *tdr* vs WT | 2,522 | 6,285 | 2,947 |
| down in both | – | 3,787 | **1,992** |

## GO annotation files

`data/GO/goid2gene_BP.txt` pairs GO BP terms with RAP locus IDs. It combines two sources:

- GO terms extracted from the RAP-DB `IRGSP-1.0_representative_annotation_2021-11-11.tsv`;
- the Oryzabase gene list (downloaded 2023-03-17).

It is restricted to BP terms using GO.db. `godb_BP.txt` gives the term names. `riceIDtable.csv` is from http://bioinformatics.fafu.edu.cn/riceidtable/.

## References

- Kashima M, Shida Y, Yamashiro T, Hirata H, Kurosaka H (2021) Intracellular and intercellular gene regulatory network inference from time-course individual RNA-Seq. *Front Bioinform* 1: 777299. https://doi.org/10.3389/fbinf.2021.777299
- Endo M, Tsuchiya T, Hamada K, Kawamura S, Yano K, Ohshima M, Higashitani A, Watanabe M, Kawagishi-Kobayashi M (2009) High temperatures cause male sterility in rice plants with transcriptional alterations during pollen development. *Plant Cell Physiol* 50: 1911–1922. https://doi.org/10.1093/pcp/pcp135
- Sun J, Nishiyama T, Shimizu K, Kadota K (2013) TCC: an R package for comparing tag count data with robust normalization strategies. *BMC Bioinformatics* 14: 219.
- Hao Y, et al. (2021) Integrated analysis of multimodal single-cell data. *Cell* 184: 3573–3587.
- Wu T, et al. (2021) clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *Innovation* 2: 100141.
- Sakai H, et al. (2013) Rice Annotation Project Database (RAP-DB). *Plant Cell Physiol* 54: e6.

## License

[to be decided, e.g. MIT]

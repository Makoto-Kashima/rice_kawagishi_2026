# OsMYB80 / TDR anther RNA-seq analysis

R script for the anther RNA-seq analysis of wild-type, *osmyb80* and *tdr* rice (*Oryza sativa* cv. Nipponbare).

> [Citation of the paper — to be added]

## Contents

```
OsMYB80_TDR_RNAseq_analysis.R   analysis script
input/                          count matrix and sample metadata (also deposited in GEO [GSExxxxxx])
```

| input file | content |
|---|---|
| `OsMYB80_TDR_anther_Salmon_NumReads.txt.gz` | Salmon-estimated read counts; rows = RAP-DB IRGSP-1.0 representative transcripts, columns = 56 samples |
| `OsMYB80_TDR_anther_sample_metadata.txt.gz` | sample information: genotype (WT n = 18, *osmyb80* n = 10, *tdr* n = 28), allele, measurements, column order |

Read processing, mapping and quantification (fastp, BWA-MEM, Salmon v1.4.0; RAP-DB IRGSP-1.0 representative transcripts, release 2022-09-01) followed Kashima et al. (2021) *Front. Bioinform.* 1:777299.

## Requirements

R 4.1 / Bioconductor 3.14: Seurat 4.2.0, SeuratObject 4.1.2, TCC 1.34.0 (edgeR 3.36.0), clusterProfiler 4.2.2, AnnotationHub 3.2.2, ggplot2, patchwork, VennDiagram.
Newer edgeR versions give different TCC results. The `bioconductor/bioconductor_docker:RELEASE_3_14` container reproduces the environment.

## Usage

```bash
Rscript OsMYB80_TDR_RNAseq_analysis.R
```

Run it from the repository root. Results are written to `results/`:
- differential expression tables (TCC; *osmyb80* vs WT and *tdr* vs WT);
- summary counts and downregulated gene sets;
- QC plots;
- `sessionInfo.txt`.

## License

MIT

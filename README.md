# multiome Pipeline

Functions to perform coembedded multiome 10X genomics analyses based in Seurat pipeline (Satija lab).
## This script include the following steps: 
- QC
- Cell filtering based on features/UMIS/perc mitochondrial DNA and also provides a way to identifycells expressing housekeeping genes as GAPHD+  to add as a filtering step if desired.
- Differential expression analyses (DEA) to identify altered genes noted as DEGs.
- Peaks reannotation using MACS3 and Differential accessibility analysis (DAA) to identify accesibility changed surrounding genes noted as DAGs.
- Peaks annotation with the closest gene.
- Venn diagram to find the intercept of DEGs and DAGs and if they follow the same regulatory sense in both snRNA-Seq and snATAC-Seq.

Developped and maintained by Mar Muniz Moreno. Last update July 2022.
This scrip have a serie of functions that will allow to perform multiome analysis in the case where a co embedded multiome is needed due to different processing times of the samples.

## For more info about the tools, tutorial and main packages used:

- https://satijalab.org/signac/articles/pbmc_multiomic.html
- https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
- https://biofam.github.io/MOFA2/tutorials.html
- https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
- https://satijalab.org/signac/articles/pbmc_multiomic.html
- https://satijalab.org/seurat/articles/multimodal_vignette.html
- https://satijalab.org/signac/articles/pbmc_vignette.html
- https://satijalab.org/signac/articles/integrate_atac.html
- https://satijalab.org/seurat/articles/atacseq_integration_vignette.html#annotate-scatac-seq-cells-via-label-transfer-1

## Libraries needed
R version 4.1.2 (2021-11-01) 

```R
library("Signac");library("Seurat");library("EnsDb.Mmusculus.v79"); library("ggplot2");
library("biovizBase"); library("hdf5r");library("dplyr");library("tidyr");library("gridExtra");
library("patchwork"); library("grid");library("gridBase"); library("cowplot");
library("BSgenome.Mmusculus.UCSC.mm10"); library("readxl");library("matrixStats");
library("VennDiagram");
```

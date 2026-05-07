# IM Cell-Fate Landscape

Analysis code and supporting data for the accepted version of:

**Single-nucleus transcriptomics resolves multiple fate dynamics between inflorescence meristem and primary stem**

Moreno-Ramirez et al., 2026. Science Advances 

An earlier preprint version is available at bioRxiv:<br>
<https://doi.org/10.1101/2024.08.06.606781>

## Overview

This repository contains Python scripts used to analyse Arabidopsis inflorescence meristem single-nucleus RNA-seq data. The workflow is built mainly around Scanpy, with Harmony used for sample integration and scFates/Palantir used for trajectory analyses.

The analyses cover:

- preprocessing of CellBender-corrected 10x count matrices;
- quality control, doublet removal, normalization, HVG selection, batch integration and Leiden clustering;
- manual cell-type annotation and marker-gene visualization;
- early primordia gene regulatory network analysis using TF-target evidence;
- cell-cycle and differentiation trajectory inference;
- gene ontology enrichment for cluster and trajectory gene sets;
- comparisons with published Arabidopsis meristem, flower and stem-cell datasets;
- validation plotting for `gh3hex` mutant measurements.

## Repository Structure

```text
.
|-- env/
|   |-- unbiased_clustering.yml
|   `-- trajectory_analysis.yml
|-- input_data/
|   |-- output_cellbender_filtered_resequenced_SITTA1.h5
|   |-- output_cellbender_filtered_resequenced_SITTC6.h5
|   `-- output_cellbender_filtered_resequenced_SITTE2.h5
|-- output_data/
|   `-- generated analysis tables and AnnData outputs
|-- figures/
|   `-- generated figures
|-- scripts/
|   |-- unbiased_clustering.py
|   |-- marker_genes_plotting.py
|   |-- GRN_EarlyPrimordia.py
|   |-- batches_analysis.py
|   |-- cell_cycle_trajectory_analysis.py
|   |-- internal_cell_types_trajectory_analysis.py
|   |-- trajectory_SC_to_EP.py
|   |-- go_analysis.py
|   |-- comparison_studies.py
|   `-- gh3hex_analysis.py

```

## Environments

Two Conda environments are provided because the clustering and trajectory workflows use different dependency sets.
For preprocessing, clustering, marker genes, GO analysis and comparison plots:
```bash
conda env create -f env/unbiased_clustering.yml
conda activate unbiased_clustering
```
For scFates/Palantir trajectory analyses:
```bash
conda env create -f env/trajectory_analysis.yml
conda activate trajectory_analysis
```

## Data input

The main inputs are the three CellBender-corrected 10x HDF5 matrices in `input_data/`. These correspond to the resequenced SAM biological samples used for the main single-nucleus analysis. Raw files are deposited in in the NCBI Sequence Read Archive under BioProject accession PRJNA1438735, with SRA experiment accessions SRX32555293, SRX32555294, and SRX32555295. 

## Main Workflow

### 1. Unbiased Clustering

Script: `scripts/unbiased_clustering.py`

This is the main preprocessing and clustering workflow. It imports the three CellBender-corrected 10x matrices, labels samples, removes predicted doublets and low-quality nuclei, removes mitochondrial/chloroplast genes, normalizes and log-transforms counts, filters outlier genes by pseudo-bulk comparison, selects HVGs, integrates samples with Harmony, performs PCA/neighbors/UMAP/Leiden clustering, assigns manual cell-type names and exports annotated AnnData objects plus marker-gene tables.

### 2. Marker-Gene Visualization

Script: `scripts/marker_genes_plotting.py`

Generates matrix plots and UMAP expression plots for manually selected marker sets, including cell-identity markers and curated hormone-related gene sets where available.

### 3. Batch and Replicate Checks

Script: `scripts/batches_analysis.py`

Plots UMAPs split by sample, compares pseudo-bulk expression between biological replicates and summarizes detected genes, total counts and cluster proportions across samples.

### 4. Early Primordia GRN

Script: `scripts/GRN_EarlyPrimordia.py`

Builds a TF-target network for the early primordia (`EP`) cluster using marker genes and DAP-seq-supported TF-target evidence.

### 5. Cell-Cycle Trajectory

Script: `scripts/cell_cycle_trajectory_analysis.py`

Uses the S-phase and G2/M-phase clusters to infer a cell-cycle trajectory with Palantir diffusion space and scFates principal tree fitting. It then identifies pseudotime-associated genes, fits expression trends and orders genes by peak expression along the trajectory.

### 6. Internal Cell-Type Differentiation Trajectories

Script: `scripts/internal_cell_types_trajectory_analysis.py`

Infers trajectories from undifferentiated cells toward internal SAM cell types, including xylem parenchyma, procambium, EP and cortex. It uses Palantir diffusion maps, MAGIC imputation and scFates tree fitting, then tests for pseudotime-associated genes along selected milestones.

### 7. Stem-Cell to Early-Primordia Trajectory

Script: `scripts/trajectory_SC_to_EP.py`

Runs a focused trajectory from undifferentiated/stem-cell identity to early primordia identity. Marker expression used for orientation includes genes such as `WUS`, `CLV3` and `AHP6`.

### 8. Gene Ontology Analysis

Script: `scripts/go_analysis.py`

Runs g:Profiler enrichment for marker genes from annotated clusters and plots top biological-process terms. This script uses the online g:Profiler service through `gprofiler-official`, so it requires internet access.

### 9. Comparison With Published Studies

Script: `scripts/comparison_studies.py`

Compares marker genes from this study with published Arabidopsis single-cell or domain-specific datasets and generates cross-study marker comparison summaries.

### 10. `gh3hex` Validation Plots

Script: `scripts/gh3hex_analysis.py`

Plots SAM size and phylotaxis/divergence-angle measurements for `gh3hex` validation experiments and runs Tukey tests with Pingouin.

## Helper Functions

`snuclei_run_def.py` contains reusable helper functions for Scrublet doublet removal, QC plotting and filtering, HVG/PCA calculation, UMAP and Leiden clustering, gene-name conversion, cell-cycle scoring, cell-domain marker summaries and related exploratory analyses.

## Suggested Run Order

For a full analysis rerun, the intended order is approximately:

1. Place the required input files in `input_data/`.
2. Run `scripts/unbiased_clustering.py` to produce the annotated AnnData object and marker tables.
3. Run `scripts/batches_analysis.py` and `scripts/marker_genes_plotting.py` for QC, replicate and marker visualizations.
4. Unzip the TF-target database if needed, then run `scripts/GRN_EarlyPrimordia.py`.
5. Run trajectory scripts in the `trajectory_analysis` environment:
   - `scripts/cell_cycle_trajectory_analysis.py`
   - `scripts/internal_cell_types_trajectory_analysis.py`
   - `scripts/trajectory_SC_to_EP.py`
6. Run `scripts/go_analysis.py` for cluster-level GO enrichment.
7. Run `scripts/comparison_studies.py` for cross-study marker comparisons.
8. Run `scripts/gh3hex_analysis.py` for validation measurements, if the validation input files are available.

## Reproducibility Notes

- The scripts were written as manuscript analysis scripts rather than as a packaged command-line pipeline.
- Most scripts now use repository-level input/output folders, but some manuscript-era sections may still require path checks before a complete rerun.
- Generated outputs can be large; keep large `.h5ad` files local or use Git LFS/data repositories when sharing them.

## Contact

Sebastian Moreno-Ramirez<br>
sebastian.moreno@slcu.cam.ac.uk
https://github.com/SebaMorenoR/

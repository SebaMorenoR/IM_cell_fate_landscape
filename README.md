# Introduction
This repository contains files associated with the manuscript:

"Single-nucleus transcriptomics reveals cellular heterogeneity and differentiation dynamics within the shoot apical meristem"

Moreno et al. (BioRxiv)
DOI: https://doi.org/10.1101/2024.08.06.606781

# Overview 
This repository contains the analysis pipeline used to process and analyse single-nucleus RNA-seq data from the shoot apical meristem (SAM).
The workflow includes:
Data preprocessing and unbiased clustering
Marker gene identification and visualization
Gene regulatory network (GRN) inference for auxin-related clusters
Cell-cycle trajectory analysis
Differentiation trajectory analysis for inner SAM cell types
Functional analysis and comparisons with published datasets

All analyses are implemented in Python using Scanpy-based workflows.

# Installation

All scripts developed in this project use Python and rely on the packages indicated in the conda environments folder.

Environment for Figures 1 to 3 in conda_environments/unbiased_clustering.yml
Example: 
* conda env create -f conda_environments/unbiased_clustering.yml conda activate unbiased_clustering

Environment for Figure 5 and 6: 
* conda env create -f conda_environments/trajectory_analysis.yml conda activate trajectory_analysis

# Scripts

To run the scripts, you first need to download the data from the University of Cambridge repository (link to be updated). Currently, all data needed to run the code is in the Data folder within this repository.

## 1.- Input Data
The analysis requires the CellBender-corrected count matrices generated from three sequencing libraries:

output_cellbender_filtered_resequenced_SITTA1.h5
output_cellbender_filtered_resequenced_SITTE2.h5
output_cellbender_filtered_resequenced_SITTC6.h5

These correspond to the three SAM biological samples used in the study.

Place these files inside the Data/ directory or modify the path variable in the scripts.

#### Output files:

* output_rank_genes_group_final_AGI.xlsx: Shows the top 150 marker genes per cluster. This is used for following analysis.
* pcHVG_data_processed.h5ad: The first generated file, which is the AnnData with HVG genes after the Leiden unbiased clustering approach.
* pcHVG_data_processed_names.h5ad: Same as before but with cluster names added.
* adata_filtered_norm_clusters.h5ad: Contains all high-quality cells and genes with PCA information necessary for trajectory analysis.


## 2.- Marker gene visualization
Script:
marker_genes_plotting.py

Generates visualizations for cluster-specific markers including:
dot plots
heatmaps
gene expression visualizations on UMAP

## 3.- Gene Regulatory network (GRN) for auxin-related cluster (Early primordia cluster) (Figure 3).
Script:
GRN_EarlyPrimordia.py

#### Inputs:
output_rank_genes_group_final_AGI.xlsx
red_oficial_pares_cruce.csv

The second file contains TF-target interaction evidence.

#### Outputs:
tf_network_pure_evidence_EP_cluster_100.xlsx
Network table compatible with Cytoscape for visualization.

The script also generates the GRN plot shown in Figure 3.

## 3.- Cell-cycle trajectory analysis (Figure 5).

### a.- Differentially expressed genes along cell-cycle trajectory.

#### Import file: 
* 'adata_filtered_norm_cluster.h5ad' (it was generated in the Unbaised clustering section. Cell-cycle clusters are isolated for trajectory analysis (S-phase and G2/M-phase clusters))

#### Output file: 
* "adata_scfates_cellcycle_fitted.h5ad" (which will be used for ploting heatmaps of differentially expressed genes and to generate "Supplementary Table", which has the heatmap DEG cell-cycle genes ordered by peak of expression )

## 4.- Differentiation trajectories from OC to primary stem (Figure 6).

### a.- Removing cell-cycle effect to increase resolution. 
* 'removing_cell_cycle_effect.py' (used to remove cell cycle effect and to assign new cluster that then will be used for trajectory )
#### Import files: 
* "adata_filtered_norm_clusters.h5ad" 
### Output files" 
* "adata_filtered_norm_clusters_without_cell_cycle.h5ad"
### b.- Trajectory analysis for internal cell-types
* "internal_cell_types_trajectory_analysis.py" function is used to obtain trajectory analysis for cortex, vasculature 1, vasculature 2, internal stem cells and early primordia clusters. 
#### Import files: 
* "adata_filtered_norm_clusters_without_cell_cycle.h5ad"
#### Output files: 

Plots observed in Figure 6 plus several excel files with the differentially expressed genes along the different trajectories. Genes are sorted by peak of expression. 

# Contact

Sebastian Moreno (sebastian.moreno@slcu.cam.ac.uk)


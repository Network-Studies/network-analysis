# network-analysis
# Protein-Protein Interaction & Pathway Enrichment Analysis

This repository contains the custom Python computational pipeline used for analyzing the biological networks of Strontium and Calcium. The pipeline constructs protein-protein interaction (PPI) networks, calculates connectivity indices, performs Monte Carlo simulations, and conducts pathway network enrichment.

### Citation & Authorship
This code accompanies the publication: **""** * **Authors:** [Insert 1st Author Name], [Insert Your Name], [Insert Other Authors]  
* **Code Development & Maintenance:** [Insert Your Name] 

## Pipeline Overview

The analysis is broken down into modular Python scripts using the ROBIO and RELIO tools:

1. **`stage1_network.py`:** Cleans input gene data and scans the HIPPIE database to extract high-confidence interacting nodes. Calculates the Gold Standard Connectivity Index (CI) and Largest Connected Component (LCC).
2. **`stage4_simulation.py`:** Performs 1,000 random Monte Carlo sampling trials on the HIPPIE network to establish a null distribution and calculates empirical p-values.
3. **`final_enrichment.py`:** Conducts hypergeometric enrichment analysis of the gene sets against the Reactome pathway database, applying FDR correction.
4. **`analyze_network.py`:** Computes advanced network metrics (Global Efficiency, Density, Centralization, Assortativity) on the significant pathways.
5. **`enrichment_plot.py`:** Generates comparative visualizations of the pathway enrichment scores between Strontium and Calcium.

## Usage

To reproduce the analysis, ensure the necessary input data files (`calcium_genes.txt`, `strontium_genes.txt`) are in the same directory. Note: The full HIPPIE database and Reactome GMT files are required to run the pipeline but are not hosted here due to size limits; they must be downloaded from their respective public repositories.

Run the scripts in the following order:
`python stage1_network.py`
`python stage4_simulation.py`
`python final_enrichment.py`
`python analyze_network.py`
`python enrichment_plot.py`


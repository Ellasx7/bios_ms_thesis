# Thesis Project: Systematic Evaluation of Spatial Transcriptomics Alignment Methods 

- This repository accompanies the thesis **“Systematic Evaluation of Spatial Transcriptomics Alignment Methods.”**  
- It provides example data, evaluation code, and a tutorial demonstrating how different alignment strategies perform when integrating spatial transcriptomics sections.

---

## Overview

This project evaluates spatial alignment methods that integrate transcriptomic and image-based features from consecutive tissue sections.  
The comparison framework supports quantitative and visual assessment of alignment accuracy based on spot-level coordinate transformations.

---

## Repository Structure

```text
project-root/
│
├── code/
│ └── eval_main.py
│
├── sample_data/
│ ├── spatial_coords/
│ │ ├── B1_to_B2_ref_gene_expression_based_coords.csv
│ │ ├── B1_to_B2_ref_image_features_based_coords.csv
│ │ ├── B1_to_B2_qry_gene_expression_based_coords.csv
│ │ └── B1_to_B2_qry_image_features_based_coords.csv
│ │
│ └── features/
│ ├── B1_to_B2_ref_gene_expression.csv
│ ├── B1_to_B2_ref_image_features.csv
│ ├── B1_to_B2_qry_gene_expression.csv
│ └── B1_to_B2_qry_image_features.csv
│
├── tutorial.ipynb
├── .gitignore
└── README.md
```

---

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/Ellasx7/ST-alignment-eval.git
cd ST-alignment-eval
```
### 2. Try the tutorial
Open the Jupyter notebook to walk through the evaluation example:
```bash
jupyter lab tutorial.ipynb
```

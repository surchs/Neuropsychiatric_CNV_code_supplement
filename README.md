# Code supplement to the Neuropsychiatric CNV FC paper
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F862615-informational
)](https://doi.org/10.1101/862615)

This repository contains the code used to process and analyse the data presented in the "Neuropsychiatric mutations delineate functional brain connectivity dimensions contributing to autism and schizophrenia" paper. 

## Abstract
16p11.2 and 22q11.2 Copy Number Variants (CNVs) confer high risk for Autism Spectrum Disorder (ASD), schizophrenia (SZ), and Attention-Deficit-Hyperactivity-Disorder (ADHD), but their impact on functional connectivity (FC) networks remains unclear. 

We analyzed resting-state functional magnetic resonance imaging data from 101 CNV carriers, 755 individuals with idiopathic ASD, SZ, or ADHD and 1,072 controls. We used CNV FC-signatures to identify major dimensions contributing to complex idiopathic conditions. 

CNVs had large mirror effects on FC at the global and regional level, and their effect-sizes were twice as large as those of idiopathic conditions. Thalamus, somatomotor, and posterior insula regions played a critical role in dysconnectivity shared across deletions, duplications, idiopathic ASD, SZ but not ADHD. Individuals with higher similarity to deletion FC-signatures exhibited worse behavioral and cognitive symptoms. 

The FC-signatures of both deletions were associated with the spatial expression pattern of genes within as well genes outside these 2 loci. This genetic redundancy may represent a factor underlying shared FC signatures between both deletions and idiopathic conditions.

![Figure 1](https://github.com/surchs/Neuropsychiatric_CNV_code_supplement/blob/master/NP_CNV_Fig1.png)

## Installation
Requirements
- [NIAK](http://niak.simexp-lab.org/build/html/index.html)
- `cnvfc` (included in this repository)
- `gene_expression` (submodule included in this repository)

To download the code of this repository, run the following command in a terminal:
```
git clone git@github.com:surchs/Neuropsychiatric_CNV_code_supplement.git --recursive
```
The `--recursive` flag will ensure that you also download the analysis code in the [`gene_expression` submodule](https://github.com/kkumar-iitkgp-livia/GeneExp_and_CNV_FCsignatures/tree/master).

The analysis scripts included in this repository make use of a number of custom python functions included in the `cnvfc` package. This package does not have to be installed but is locally referenced.

## How to use the repository
The notebooks in `Notebooks/16p_FC_profile.ipynb` and `Notebooks/22q_FC_profile.ipynb` illustrate the identified FC signatures for the 16p11.2 and 22q11.2 deletion carriers respectively. These findings are based on the analysis scripts in the `Scripts` folder that have been run in the following order:

1.   `Scripts/preprocess_data.m` is an example `NIAK` preprocessing script to preprocess the raw anatomical and functional data.
2.   `Scripts/generate_connectomes.m` is an example `NIAK` analysis script to compute the seed-based functional connectome for the MIST_64 atlas. This analysis step also implements the regression of noise confounds from the preprocessed functional time series data.
3.   `Scripts/recast_vectorized_connectome_to_matrix.py` is a helper script to re-organize the vectorized connectome from the Matlab style column major ordering to the numpy style row major ordering.
4.   `Scripts/FC_case_control_contrast.py` is an analysis script to compute the case-control FC profiles reported in the paper (including the 16p11.2 and 22q11.2 CNV FC profiles)
5.   `Scripts/CNV_FC_profile_enrichment.py` is an analysis script to compute the similarity between the FC connectomes of individuals in the idiopathic neuropsychiatric samples and the identified CNV FC profiles.
6.   `Scripts/null_model.py` is a helper script to compute a distribution of randomly permuted FC profiles to compute exact p-values against.
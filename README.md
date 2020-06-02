# Code supplement to the Neuropsychiatric CNV FC paper
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![DOI](https://github.com/surchs/Neuropsychiatric_CNV_code_supplement/blob/master/cnv_paper_badge.svg
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

The analysis scripts included in this repository make use of a number of custom python functions included in the `cnvfc` package. This package does not have to be installed but is locally referenced.

# Simulated GWAS & Statistical Power Analysis

A Python notebook implementation that simulates a quantitative Genome-Wide Association Study (GWAS), runs single-SNP linear association testing, checks genomic inflation, and evaluates statistical power across sample sizes.

## Overview

The notebook models genotype and continuous phenotype data synthetically to verify how association statistics, multiple-testing correction, and sample size affect causal variant discovery. It includes:

* **Synthetic GWAS simulation (`simulate_gwas`)** generating minor allele frequencies (`maf`), diploid genotypes (`0/1/2`), standardized genotype dosages, sparse causal SNP effects (`beta`), genetic liability scaled to target heritability (`h2`), and continuous traits (`y`)
* **Fast vectorized association testing (`association`)** computing Pearson correlations (`r`), `t`-statistics, and two-sided `p`-values against standardized phenotypes
* **Diagnostic & discovery metrics** reporting the genomic inflation factor ($\lambda_{\text{GC}}$), Bonferroni significance threshold ($\alpha = 0.05 / M$), total significant hits, and true causal SNPs recovered
* **Visualizations** producing a chromosome-styled Manhattan plot (`-log10(p)`), a QQ plot against expected uniform quantiles, and a log-scale sample-size power curve



## Default Simulation & Results

With the default parameters (`N = 2000` individuals, `M = 5000` SNPs across 5 simulated chromosomes, `n_causal = 5`, `h2 = 0.25`, `seed = 42`), the notebook outputs:

* **Genomic inflation ($\lambda_{\text{GC}}$):** `0.97`
* **Bonferroni threshold:** `p < 1.0e-05`
* **Significant SNPs:** `5`
* **Causal SNPs recovered:** `5 of 5` (SNPs `1048`, `1108`, `1244`, `2769`, and `3455`)



### Sample Size vs. Power

Running 10 replicates per sample size (`N` from `200` to `5000` with `seed = 1000 + rep` and default `M = 5000`, `n_causal = 5`, `h2 = 0.25`) yields:

| Sample Size (`N`) | Mean Power (Causal SNPs Detected) | Mean False Positives |
| --- | --- | --- |
| 200 | 0.16 | 0.1 |
| 500 | 0.62 | 0.0 |
| 1000 | 0.98 | 0.0 |
| 2000 | 1.00 | 0.0 |
| 5000 | 1.00 | 0.2 |

## Requirements

* Python 3
* numpy
* pandas
* matplotlib
* scipy


## Usage

1. Open the notebook in Google Colab or Jupyter Notebook.
2. Run the simulation and association cell to generate the SNP results table and summary metrics.
3. Run the plotting cells to inspect the Manhattan plot, QQ plot, and power-versus-sample-size curve.

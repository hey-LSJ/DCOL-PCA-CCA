# Fast and Tuning-free Nonlinear Embedding and Integration for Omics Data

## Background

The rapid advancement of single-cell technologies has led to an explosion of diverse omics data, significantly enhancing our understanding of cell populations, disease states, and developmental lineages. However, traditional dimension reduction (DR) techniques, particularly linear methods like PCA, often struggle to capture the complex associations inherent in this data. While nonlinear techniques show promise, they frequently encounter scalability issues and may not effectively focus on informative embeddings.

## Method

Here, we introduce **DCOL-PCA** and **DCOL-CCA**, our innovative methods for dimension reduction and integration of single- and multi-omics data. By leveraging the **Dissimilarity based on Conditional Ordered List (DCOL)** correlation, these methods quantify nonlinear relationships between variables, effectively addressing the challenges posed by high-dimensional, noisy, and sparse single-cell omics data.

## Features
- **DCOL-PCA**: A method for dimension reduction based on DCOL correlation. It takes an input omics matrix and returns a lower-dimensional representation that preserves essential biological features. By effectively capturing nonlinear relationships, DCOL-PCA is robust against high-dimensional, noisy, and sparse data commonly encountered in single-cell omics.

- **DCOL-CCA**: An innovative method for integrating multi-omics data based on DCOL correlation. This method enables researchers to align and correlate datasets from different omics layers, effectively capturing the complex interactions between them. DCOL-CCA takes two input omics matrices, X (N x p) and Y (N x q), and returns lower-dimensional representations for each layer, respectively.

## Tutorial
We provide examples on how to use **DCOL-PCA** and **DCOL-CCA**. Their usage is straightforward. Please refer to the [R Markdown tutorial](https://github.com/hey-LSJ/DCOL-PCA-CCA/blob/main/tutorial/tutorial.Rmd).




For further details, please refer to the full paper or contact the authors for inquiries.

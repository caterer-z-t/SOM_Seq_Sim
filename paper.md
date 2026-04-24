---
title: 'SOM-Seq: A Python Toolbox for Single-Cell Sequencing Simulation and Self-Organizing Map Analysis'
tags:
  - Python
  - bioinformatics
  - single-cell sequencing
  - self-organizing maps
  - dimensionality reduction
  - clustering
  - simulation
authors:
  - name: Zachary Caterer
    orcid: 0000-0001-9019-0730
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2"
  - name: Madeline Pernat
    orcid: 0000-0003-2814-3428
    equal-contrib: true
    corresponding: true
    affiliation: 3
  - name: Victoria Hurd
    orcid: 0000-0002-5548-6883
    equal-contrib: true
    corresponding: true
    affiliation: 4
affiliations:
  - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO, USA
    index: 1
  - name: Department of Biomedical Informatics, University of Colorado Anschutz Medical Campus, Aurora, CO, USA
    index: 2
  - name: Department of Civil, Environmental, and Architectural Engineering, University of Colorado Boulder, Boulder, CO, USA
    index: 3
  - name: Ann and H.J. Smead Department of Aerospace Engineering Sciences, University of Colorado Boulder, Boulder, CO, USA
    index: 4
date: 24 April 2026
bibliography: paper.bib
---

# Summary

SOM-Seq is an open-source Python toolbox that integrates two complementary workflows for single-cell genomics research: synthetic sequencing data generation (`Seq_Sim`) and Self-Organizing Map (SOM) based clustering and visualization (`SOM`). The `Seq_Sim` module generates realistic pseudo-bulk single-cell datasets with configurable cell-type compositions, batch effects, disease states, and differential expression patterns, adapted from simulation approaches developed in the Zhang Lab [@inamo2024scorpio]]. The `SOM` module provides a high-level Python class built on MiniSom [@vettigli2018minisom] that handles data scaling, automated hyperparameter tuning, topographic quality metrics, and publication-quality visualizations including component planes and categorical overlays. Both modules expose command-line interfaces (CLIs), making them composable within broader bioinformatics pipelines or usable as standalone tools.

# Statement of Need

High-dimensional single-cell RNA sequencing (scRNA-seq) data requires dimensionality reduction and clustering to reveal biologically meaningful structure [@luecken2019]. Established tools such as Seurat [@hao2021] and Scanpy [@wolf2018] are now standard in single-cell analysis, typically pairing graph-based community detection with t-SNE [@van_der_maaten2008] or UMAP [@mcinnes2018] for visualization. While powerful, these methods embed data into a continuous low-dimensional space that does not explicitly preserve the topological distances between clusters, making it difficult to reason about the relative proximity of cell populations.

Self-Organizing Maps (SOMs), introduced by Kohonen [@kohonen1990], address this limitation by producing a discrete two-dimensional grid of neurons in which neighboring neurons represent similar regions of the input feature space, explicitly encoding topological structure. Despite their interpretive advantages, SOMs remain underutilized in the single-cell community. A key barrier is the absence of a well-tested, end-to-end Python package that combines SOM fitting with single-cell-style simulation, reducing the friction required to benchmark SOM-based clustering against other methods on controlled synthetic data.

SOM-Seq addresses this gap in two ways. First, `Seq_Sim` generates synthetic datasets that statistically mirror real single-cell data—including heterogeneous cell-type proportions, batch variability, and disease-associated fold changes—providing researchers with a reproducible, ground-truth-labeled environment for method comparison without requiring access to patient data. Second, the `SOM` module wraps the complete SOM workflow (scaling, training, metric evaluation, and visualization) into a clean Python API and CLI, lowering the expertise required to apply SOM-based analysis to tabular omics data.

Together, these modules enable researchers to simulate a dataset with known structure, fit a SOM, and immediately evaluate clustering quality using Percent Variance Explained (PVE) and topographic error—a capability not available in existing single-cell analysis frameworks.

# Design and Implementation

## Sequence Simulation (`Seq_Sim`)

The `Seq_Sim` module generates synthetic single-cell datasets by constructing a subject-level metadata table (age, sex, disease status, batch) and a cell-type composition matrix. Cell counts for major and rare cell populations are drawn from uniform distributions parameterized by user-supplied standard deviations and relative abundances. Disease-associated differential abundance is introduced by adding or removing cells of specified types in proportion to a configurable fold-change parameter. Pseudo-feature expression matrices are then generated per cell by combining cluster-specific signal, disease variance, and individual-level variance, with additive Gaussian noise controlled by a cluster ratio parameter. All random operations accept a seed argument to guarantee reproducibility.

Key configurable parameters exposed via `config.yml` or CLI include:

- `num_samples`: number of subjects to simulate
- `fold_change`: magnitude of disease-associated differential abundance
- `n_major_cell_types`, `n_minor_cell_types`: cell-type composition
- `n_features`: number of pseudo-expression features per cell
- `n_batches`, `prop_disease`, `prop_sex`: study design parameters

## Self-Organizing Maps (`SOM`)

The `SOM` module wraps MiniSom [@vettigli2018minisom] into a `SOM` Python class that manages the full analysis workflow:

1. **Input validation**: type and range checks on all constructor arguments.
2. **Scaling**: z-score or min-max normalization with corresponding inverse transforms stored for weight unscaling.
3. **Training**: configurable grid dimensions, rectangular or hexagonal topology, and Gaussian or bubble neighborhood functions.
4. **Metrics**: PVE measures how much of the input variance is captured by the neuron weight vectors; topographic error [@kohonen1990] quantifies how often a data point's two best-matching units are non-adjacent on the grid.
5. **Hyperparameter tuning**: when multiple values are supplied for any hyperparameter, the CLI evaluates all combinations and selects the configuration maximizing `PVE − 100 × topographic_error`.
6. **Visualization**: component plane plots (one heatmap per input feature) and categorical overlay plots (one heatmap per categorical variable), saved as publication-quality PNG files.

# Usage

## Generating Sequencing Data

```bash
python Seq_Sim/seq_sim.py \
    --num_samples 30 \
    --fold_change 0.5 \
    --config_file Seq_Sim/config.yml
```

## Fitting a SOM

```bash
python SOM/som.py \
    -t data/seq_sim_training_data.csv \
    -c data/seq_sim_categorical_data.csv \
    -o output/ \
    -s zscore \
    -x 5 -y 4 \
    -p hexagonal \
    -n gaussian \
    -e 100
```

## Python API

```python
import pandas as pd
from SOM.utils.som_utils import SOM

train = pd.read_csv("data/seq_sim_training_data.csv")
meta  = pd.read_csv("data/seq_sim_categorical_data.csv")

som = SOM(
    train_dat=train,
    other_dat=meta,
    scale_method="zscore",
    x_dim=5,
    y_dim=4,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=100,
)
som.train_map()
print(f"PVE: {som.calculate_percent_variance_explained():.1f}%")
print(f"Topographic error: {som.calculate_topographic_error():.3f}")
som.plot_component_planes(output_dir="output/")
som.plot_categorical_data(output_dir="output/")
```

# Testing and Documentation

SOM-Seq ships with a `pytest` test suite covering both modules, including input validation, scaling round-trips, training, metric calculations, and plot generation. Continuous integration via GitHub Actions runs the full suite on each push. API documentation is hosted on Read the Docs.

# Acknowledgements

<!-- TODO: Add acknowledgements (funding, advisors, etc.) -->

# References

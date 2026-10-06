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
    affiliation: "1, 2, 3"
  - name: Madeline Pernat
    orcid: 0000-0003-2814-3428
    equal-contrib: true
    corresponding: true
    affiliation: 4
  - name: Victoria Hurd
    orcid: 0000-0002-5548-6883
    equal-contrib: true
    corresponding: true
    affiliation: 5
affiliations:
  - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO, USA
    index: 1
  - name: Department of Biomedical Informatics, University of Colorado Anschutz Medical Campus, Aurora, CO, USA
    index: 2
  - name: Biofrontiers Institute Interdisciplinary Biology PhD Program, University of Colorado Boulder, Boulder, CO, USA
    index: 3
  - name: Department of Civil, Environmental, and Architectural Engineering, University of Colorado Boulder, Boulder, CO, USA
    index: 4
  - name: Department of Aerospace Engineering Sciences, University of Colorado Boulder, Boulder, CO, USA
    index: 5
date: 24 April 2026
bibliography: paper.bib
---

# Summary

Single-cell sequencing technologies profile gene expression in individual cells, enabling researchers to identify distinct cell populations and understand how they change in disease. Analyzing these datasets requires dimensionality reduction and clustering methods to reveal biologically meaningful structure. Self-Organizing Maps (SOMs) are a class of unsupervised neural network that project high-dimensional data onto a two-dimensional grid of neurons, where neighboring neurons represent similar regions of the input space, explicitly encoding topological relationships between clusters that standard methods such as UMAP or t-SNE do not preserve. Despite this interpretive advantage, SOMs remain underutilized in the single-cell community, in part because no existing tool combines SOM-based analysis with single-cell-style simulation in a single package. Researchers wishing to benchmark SOM-based clustering must currently assemble a simulation tool and a separate SOM library themselves.

SOM-Seq is an open-source Python toolbox that addresses this gap by integrating two complementary workflows: simulated single-cell dataset generation (`Seq_Sim`) and SOM-based clustering and visualization (`SOM`). The `Seq_Sim` module generates realistic simulated single-cell datasets — producing per-cell feature matrices with configurable cell-type compositions, batch effects, disease states, and differential expression patterns — providing a reproducible, ground-truth-labeled environment for benchmarking clustering methods without requiring access to patient data. This simulation approach is adapted from methods developed in the Zhang Lab [@inamo2024scorpio]. The `SOM` module provides a high-level Python class built on MiniSom [@vettigli2018minisom] that handles data scaling, automated hyperparameter tuning, topographic quality metrics, and publication-quality visualizations including component planes and categorical overlays. Both modules expose command-line interfaces (CLIs), making them composable within broader bioinformatics pipelines or usable as standalone tools.

# Statement of Need

Self-Organizing Maps (SOMs), introduced by Kohonen [@kohonen1990], produce a discrete two-dimensional grid of neurons in which neighboring neurons represent similar regions of the input feature space, explicitly encoding topological structure between clusters. Despite this interpretive advantage, SOMs remain underutilized in the single-cell community. A key barrier is the absence of a well-tested, end-to-end Python package that combines SOM fitting with single-cell-style simulation, reducing the friction required to benchmark SOM-based clustering against other methods on controlled synthetic data.

SOM-Seq addresses this gap in two ways. First, `Seq_Sim` generates synthetic datasets that statistically mirror real single-cell data—including heterogeneous cell-type proportions, batch variability, and disease-associated fold changes—providing researchers with a reproducible, ground-truth-labeled environment for method comparison without requiring access to patient data. Second, the `SOM` module wraps the complete SOM workflow (scaling, training, metric evaluation, and visualization) into a clean Python API and CLI, lowering the expertise required to apply SOM-based analysis to tabular omics data.

Together, these modules enable researchers to simulate a dataset with known structure, fit a SOM, and immediately evaluate clustering quality using Percent Variance Explained (PVE) and topographic error—a capability not available in existing single-cell analysis frameworks.

# State of the Field

High-dimensional single-cell RNA sequencing (scRNA-seq) data requires dimensionality reduction and clustering to reveal biologically meaningful structure [@luecken2019]. Established tools such as Seurat [@hao2021] and Scanpy [@wolf2018] are the standard for single-cell analysis, typically pairing graph-based community detection with t-SNE [@van_der_maaten2008] or UMAP [@mcinnes2018] for visualization. While powerful, these methods embed data into a continuous low-dimensional space that does not explicitly preserve the topological distances between clusters, making it difficult to reason about the relative proximity (i.e., similarity) of cell populations.

At the algorithmic level, general-purpose SOM libraries such as MiniSom [@vettigli2018minisom] expose raw training routines but provide no single-cell-specific workflow: users must independently implement data scaling, select grid dimensions and neighborhood parameters, compute quality metrics such as topographic error and percent variance explained (PVE), and build visualizations suited to omics data. Similarly, single-cell simulators such as the SCORPIO framework [@inamo2024scorpio], on which `Seq_Sim` is based, are not paired with an SOM analysis workflow. While a motivated user could combine these tools manually, doing so requires non-trivial implementation effort: writing a hyperparameter search over grid dimensions and neighborhood functions, implementing the PVE and topographic error calculations, and building component plane and categorical overlay visualizations — none of which MiniSom or SCORPIO provide. SOM-Seq packages all of these steps behind a single, consistent API and CLI, reducing the expertise and implementation effort required to apply SOM-based analysis to single-cell data. To our knowledge, no existing Python package integrates this complete SOM workflow, scaling, tuning, quality metrics, and omics-specific visualization, with single-cell-style data simulation in one tool.

# Software Design

SOM-Seq's design reflects two deliberate trade-offs. First, rather than requiring real patient-derived single-cell data, `Seq_Sim` generates statistically realistic synthetic data with known ground-truth cell-type labels; this sacrifices biological realism for reproducibility and removes the privacy and data-access barriers that complicate benchmarking clustering methods on human subject data. Second, rather than reimplementing SOM training, the `SOM` module wraps the existing, well-tested MiniSom library [@vettigli2018minisom] and adds the scaling, tuning, metric, and visualization steps that a complete single-cell workflow requires; this keeps the core training algorithm maintained upstream while concentrating SOM-Seq's engineering effort on the parts of the workflow specific to omics data analysis.

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
python Seq_Sim/seq_sim.py --num_samples 30 --fold_change 0.5 --config_file Seq_Sim/config.yml
```

## Fitting a SOM

```bash
python SOM/som.py -t data/sim_data_pseudo_feature_num_samples_30_fc_0.5.csv -c data/sim_data_latent_data_num_samples_30_fc_0.5.csv -o output/ -s zscore -x 5 -y 4 -p hexagonal -n gaussian -e 100
```

## Python API

```python
import pandas as pd
from SOM.utils.som_utils import SOM

train = pd.read_csv("data/sim_data_pseudo_feature_num_samples_30_fc_0.5.csv")
meta  = pd.read_csv("data/sim_data_latent_data_num_samples_30_fc_0.5.csv")

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

# Research Impact Statement

Beyond its origin as a course project, the `SOM` module is actively being applied in ongoing biomedical informatics research as a topology-preserving alternative to graph-based and explainable ML clustering pipelines. Specifically, it is being evaluated on publicly available single-cell datasets — COVID-19 PBMC multi-omics data [@stephenson2021] and ulcerative colitis tissue transcriptomics [@smillie2019intra] — as a complementary benchmarking tool alongside CellPhenoX [@young2025cellphenox], a published explainable machine learning method for single-cell clinical phenotyping. The research question being addressed is whether topology-preserving SOM-based clustering identifies biologically coherent cell populations consistent with those identified by CellPhenoX on the same datasets. Reproducible example analyses are provided in the repository: `Seq_Sim/walkthrough/seq_sim_workflow.ipynb` demonstrates the full data simulation workflow, and `SOM/examples/` contains three worked examples applying the `SOM` module to the Iris dataset (`som_example_iris.ipynb`), the Titanic dataset (`som_example_titanic.ipynb`), and simulated single-cell data (`som_example_seq.ipynb`). Community-readiness is further supported by a `pytest` test suite exercising both modules, continuous integration on every push, published API documentation, an OSI-approved open-source license, and packaging for installation via `PyPI`.

# Acknowledgements

This work originated as part of the University of Colorado Boulder course **CSCI 6118: Software Engineering for Scientists**, created by **Dr. Ryan Layer** and taught by **Dr. Erik Johnson**. We thank both Dr. Layer and Dr. Johnson for their guidance and for creating the collaborative environment that made this project possible. This work was also supported in part by the U.S. National Science Foundation under **Award No. 2022138**.

# AI usage disclosure

We used Anthropic's Claude (Sonnet 4.6) to assist with refining documentation and addressing reviewer comments, including adding clarity, removing redundant information, and drafting revised paper sections. We also used Anthropic's Claude (Sonnet 5, via Claude Code) during pre-submission preparation to identify and fix two command-line interface bugs (an argument-parsing crash in the `SOM` CLI and an invalid function call in the `Seq_Sim` CLI), to add continuous integration workflows for package build validation and PyPI publishing, and to correct inconsistencies in the LICENSE and README. During revision, Claude (Sonnet 4.6) assisted in identifying and fixing four additional bugs: an incorrect topographic error calculation for hexagonal grids, a broken fold-change effect in the simulation module, missing support for optional metadata in the CLI and API, and a non-functional `-m` plot suppression flag. All AI-assisted content and code changes were reviewed by the authors, who validated correctness and made all core design decisions.

# References

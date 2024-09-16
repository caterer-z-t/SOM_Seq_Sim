# SOM_Seq_Sim

Implementation for Self-Organizing Maps on Simulated Genetic Sequencing Datasets

## Title
Simulation of Genetic Sequencing Data using Statistically Similar Distributions from Original Datasets

## Releases & Builds

![GitHub release (with filter)](https://img.shields.io/github/v/release/caterer-z-t/SOM_Seq_Sim)
![GitHub Release Date - Published_At](https://img.shields.io/github/release-date/caterer-z-t/SOM_Seq_Sim)

## Course:
- [CSCI 6118](https://github.com/swe4s): Software Engineering for Scientists at the [University of Colorado Boulder](https://www.colorado.edu/)

## Project Overview
This project aims to develop a method for simulating genetic sequencing data using Self-Organizing Maps (SOMs). We leverage SOMs to capture and model the underlying statistical distributions of genetic data and use them to generate new datasets that are statistically similar to the original genetic sequences. This approach could be highly useful in areas like population genetics, personalized medicine, and genomic data analysis, where having access to large, synthetic, yet realistic datasets is crucial.

### Key Objectives:
- Self-Organizing Maps (SOMs): Implement SOMs to learn the structure of genetic sequencing data.
- Simulation: Develop a method to generate new genetic sequences based on statistically similar distributions learned from the SOMs.
- Evaluation: Measure the statistical similarity between the simulated and original genetic sequencing data using metrics such as Kullback–Leibler divergence and others.

## Table of Contents
1. [Project Setup](#project-setup)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Project Structure](#project-structure)
5. [Methodology](#methodology)
6. [Technologies Used](#technologies-used)
7. [Results and Evaluation](#results-and-evaluation)
8. [Future Work](#future-work)
9. [Contributions](#contributions)
10. [License](#license)

## Project Setup
### Requirements
The project requires the following software and libraries: 
- see the `env.yml` file.

## Installation
1. Clone the repository:

``` bash
git clone https://github.com/caterer-z-t/SOM_Seq_Sim.git
cd  SOM_Seq_Sim
```

2. Create a virtual environment (optional but recommended):

If you're using `conda` or `mamba`, you can create the environment directly from the `env.yml` file:
``` bash 
conda env create -f env.yml
```
Alternatively if you are using `mamba` for faster environment solving:
``` bash 
mamba env create -f env.yml
```

3. Activate environment:
``` bash
conda activate som-sim-env
```

4. Pull example data from github into Seq_Sim Directory:
``` bash
cd Seq_Sim
git clone https://github.com/borenstein-lab/microbiome-metabolome-curated-data.git
```

5. Run the SOM genetic sequencing simulation:
``` bash 
python som_genetic_simulation.py --input data/original_genetic_data.csv
```

## Usage
Command-line Arguments

    --input : Path to the original genetic sequencing dataset (CSV file).
    --output : Path to store the simulated data (default: output/simulated_genetic_data.csv).
    --som-size : Size of the SOM grid (e.g., 20x20).

### Example:
``` bash 
python som_genetic_simulation.py --input data/original_genetic_data.csv --output data/simulated_data.csv --som-size 20x20
```

## Project Structure
``` bash 
SOM_Seq_Sim/
├── docs/
│   ├── img/                     # Contains images used in Read the Docs
│   ├── pipelines/               # Additional MD files for pipelines and packages
│   ├── index.md                 # Initial page for Read the Docs
│   ├── requirements.txt         # Dependencies for Read the Docs 
├── Seq_Sim/
│   ├── dev_sequencing_code.py   # Main script to run the simulation (refer to README.md for usage)
│   ├── utils.py                # Helper functions (data processing, evaluation)
│   ├── unimodal_example.py     # Example use for unimodal distribution
│   ├── bimodal_example.py      # Example use for bimodal distribution
│   └── trimodal_example.py     # Example use for trimodal distribution
├── test/
│   └── seq_test.py             # Testing script for utils functions in Seq_Sim
├── .gitignore                  # Git ignore file to exclude files and directories from Git
├── .readthedocs.yaml           # YAML file specifying dependencies for Read the Docs
├── changelog.md                # MD file documenting changes made in the development process
├── env.yml                     # YAML file specifying dependencies for the entire package
├── README.md                   # README file
└── LICENSE                     # License for the project
```

## Methodology
### Self-Organizing Maps (SOMs)
Self-Organizing Maps are a type of artificial neural network trained using unsupervised learning to produce a lower-dimensional (typically 2D) representation of high-dimensional data. In this project, SOMs are used to learn and capture the structure of genetic sequencing data, which is then used to generate statistically similar data points.

### Simulation Approach
1. Data Preprocessing: The original genetic data is preprocessed (e.g., normalization, feature selection) before being fed into the SOM.
Training the SOM: The SOM learns the topological relationships between different genetic sequences.
2. Data Generation: New synthetic data points are generated based on the learned distribution.
3. Evaluation: The synthetic data is compared with the original dataset using statistical similarity metrics.

## Technologies Used
- Python: Programming language for implementing the algorithm.
- NumPy/Pandas: Libraries for data processing.
- SOMPY: A Python library for creating and training Self-Organizing Maps.
- Scikit-learn: Used for auxiliary machine learning tasks and model evaluation.
- Matplotlib: For visualizing the results of the SOM and data simulation.

## Results and Evaluation
We evaluate the synthetic genetic sequences by comparing them to the original dataset using various statistical measures, including:

- Kullback-Leibler Divergence: Measures the difference between probability distributions.
- Mean Squared Error (MSE): Compares the original and simulated genetic data features.

Preliminary results suggest that the simulated data retains significant statistical similarity to the original genetic dataset, indicating that SOMs can effectively learn the underlying distribution.

## Future Work
- Improved Feature Representation: Experiment with different feature representations of genetic data to improve the fidelity of the simulation.
- Larger Datasets: Apply the method to larger and more diverse genetic sequencing datasets.
- Evaluation Metrics: Explore additional metrics to assess the quality of the generated data.

## Contributors

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/madrpernat">
        <img src="https://github.com/madrpernat.png" width="200" />
      </a>
      <br />
      <b>Maddie</b>
    </td>
    <td align="center">
      <a href="https://github.com/pnnoc">
        <img src="https://github.com/pnnoc.png" width="200" />
      </a>
      <br />
      <b>Nattapat</b>
    </td>
    <td align="center">
      <a href="https://github.com/victoria-hurd">
        <img src="https://github.com/victoria-hurd.png" width="200" />
      </a>
      <br />
      <b>Viki</b>
    </td>
    <td align="center">
      <a href="https://github.com/caterer-z-t">
        <img src="https://github.com/caterer-z-t.png" width="200" />
      </a>
      <br />
      <b>Zac</b>
    </td>
  </tr>
</table>

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## Acknowledgements

We would like to thank the contributors and open-source community for their valuable contributions to this project.

<!-- 

This image will not work as long as the repo is private, commenting out for now

<a href="https://github.com/caterer-z-t/SOM_Sim_Seq/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=caterer-z-t/SOM_Sim_Seq" />
</a>
-->

## Changelog
Please refer to the `changlog.md` file for any updates and changes

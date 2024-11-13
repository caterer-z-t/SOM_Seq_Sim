# SOM_Seq_Sim

<img src="assets/som_package_logo.png" alt="SOM SIM SEQ Logo" width="25%">

Implementation for Self-Organizing Maps on Simulated Genetic Sequencing Datasets

## Title
Simulation of Genetic Sequencing Data using Statistically Similar Distributions from Original Datasets

## Releases & Builds

![GitHub release (with filter)](https://img.shields.io/github/v/release/caterer-z-t/SOM_Seq_Sim)
![GitHub Release Date - Published_At](https://img.shields.io/github/release-date/caterer-z-t/SOM_Seq_Sim)
![Tests](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/test.yml?label=tests&logo=github)

## Course:
- [CSCI 6118](https://github.com/swe4s): Software Engineering for Scientists at the [University of Colorado Boulder](https://www.colorado.edu/)

## Project Overview
This project aims to develop a method for simulating genetic sequencing data using Self-Organizing Maps (SOMs). We leverage SOMs to capture and model the underlying statistical distributions of genetic data and use them to generate new datasets that are statistically similar to the original genetic sequences. This approach could be highly useful in areas like population genetics, personalized medicine, and genomic data analysis, where having access to large, synthetic, yet realistic datasets is crucial.

## Table of Contents
1. [Project Setup](#project-setup)
2. [Installation](#installation)
3. [Project Structure](#project-structure)
4. [Contributions](#contributions)
5. [License](#license)

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

If you're using `conda`, `mamba`, or `micromamba`, you can create the environment directly from the `env.yml` file:
``` bash 
conda env create -f env.yml
```
If you are using `mamba` for environment solving:
``` bash 
mamba env create -f env.yml
```
If you are using `micromamba` for environment solving:
``` bash
micromamba env create -f env.yml
```

3. Activate environment using conda (please change {conda} to be whatever used to create the environemnt):
``` bash
conda activate som-sim-env
```

4. For generating simulated sequencing data please see `Seq_Sim/README.md` for more information. 


## Project Structure
``` bash 
SOM_Seq_Sim/
├── docs/
│   ├── img/                     # Contains images used in Read the Docs
│   ├── pipelines/               # Additional MD files for pipelines and packages
│   ├── index.md                 # Initial page for Read the Docs
│   └── requirements.txt         # Dependencies for Read the Docs 
├── SOM/
├── Seq_Sim/
│   ├── data/                   # Subsequent datasets generated from running the sequencing code
│   ├── figures/                # Figures generated from running the sequencing code
│   ├── utils/                  # folder containg additional utility functions
│   │   ├── __init__.py
│   │   └── function.R          # contains all necessary functions for simulation.R
│   ├── __init__.py
│   ├── config.yml              # configuration file use to specify all components in simulation.R
│   ├── README.md               # information specific to running this simulation sequencing generation code
│   ├── simulation.R            # main function which runs the seqeucning generation code
│   └── simulation.sh           # shell script for running sequencing generating code multiple times, see Seq_Sim/README.md for more information
├── test/
├── .gitignore                  # Git ignore file to exclude files and directories from Git
├── .readthedocs.yaml           # YAML file specifying dependencies for Read the Docs
├── changelog.md                # MD file documenting changes made in the development process
├── env.yml                     # YAML file specifying dependencies for the entire package
├── README.md                   # README file
└── LICENSE                     # License for the project
```

## Contributors

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/madrpernat">
        <img src="https://github.com/madrpernat.png" width="200" />
      </a>
      <br />
      <b>Maddy</b>
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
      <b>Vicki</b>
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

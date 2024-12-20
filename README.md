<img src="assets/som_package_logo.png" alt="SOM SIM SEQ Logo" width="25%">

# **Self-Organizing Maps to Explore Genetic Sequencing Simulation Data (SOM-Seq)**

[![GitHub release](https://img.shields.io/github/v/release/caterer-z-t/SOM_Seq_Sim)](https://github.com/caterer-z-t/SOM_Seq_Sim/releases) 
[![GitHub Release Date](https://img.shields.io/github/release-date/caterer-z-t/SOM_Seq_Sim)](https://github.com/caterer-z-t/SOM_Seq_Sim/releases) 
[![Tests](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/test.yml?label=tests&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/actions)
[![Build Status](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/build.yml?label=build&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/actions)
[![License](https://img.shields.io/github/license/caterer-z-t/SOM_Seq_Sim?label=License)](https://github.com/caterer-z-t/SOM_Seq_Sim/blob/main/LICENSE)
[![Repo Size](https://img.shields.io/github/repo-size/caterer-z-t/SOM_Seq_Sim?label=Repo%20Size)](https://github.com/caterer-z-t/SOM_Seq_Sim)
[![codecov](https://codecov.io/gh/<your-username>/<your-repo>/branch/main/graph/badge.svg)](https://codecov.io/gh/caterer-z-t/SOM_Seq_Sim/tree/zac)
[![Contributors](https://img.shields.io/github/contributors/caterer-z-t/SOM_Seq_Sim?label=Contributors)](https://github.com/caterer-z-t/SOM_Seq_Sim/graphs/contributors)
[![Pull Requests](https://img.shields.io/github/issues-pr/caterer-z-t/SOM_Seq_Sim?label=Pull%20Requests)](https://github.com/caterer-z-t/SOM_Seq_Sim/pulls)
[![Discussions](https://img.shields.io/github/discussions/caterer-z-t/SOM_Seq_Sim?label=Discussions&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/discussions)
[![Docs](https://img.shields.io/badge/docs-online-blue?logo=readthedocs)](https://som-seq-sim.readthedocs.io/en/latest/index.html)
[![Stars](https://img.shields.io/github/stars/caterer-z-t/SOM_Seq_Sim?label=Stars&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/stargazers)

---

## **Course:**

[CSCI 6118](https://github.com/swe4s): Software Engineering for Scientists at the [University of Colorado Boulder](https://www.colorado.edu/)

---

## **Project Overview**

Welcome to **SOM-Seq** 🎉, a versatile toolbox designed to tackle two distinct but complementary tasks: **sequence data generation/simulation** and the use of **Self-Organizing Maps (SOMs)** for clustering and visualizing high-dimensional data. By combining these efforts, this project explores the applicability and value of SOMs as a powerful tool for analyzing and uncovering patterns within complex sequencing data. Additionally, each part of this project can be used independently, catering to those who are specifically interested in either generating sequence data or leveraging SOMs for data analysis. 🚀

---

## **Table of Contents**

- [Usage](#usage)
  - [Project Setup](#project-setup)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Generating Sequencing Data](#generating-sequencing-data)
  - [Running SOMs](#running-soms)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Code of Conduct](#code-of-conduct)
- [Testing](#testing)
- [Contributors](#contributors)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Project Structure](#project-structure)

---

## **Usage**

### **Project Setup**

#### **Requirements**

This project depends on several libraries and software packages. You can find the list in the `env.yml` file, which makes setting up the environment a breeze.

> Key libraries in the environment:
>   - `numpy`  
>   - `pandas`  
>   - `matplotlib`  
>   - `seaborn`  
>   - `scikit-learn`  
>   - `scipy`  
>   - `statsmodels`  
>   - `simpy`  
>   - `jupyter`  
>   - `jupyterlab`  
>   - `black`
>   - `mkdocs` 
>   - `pytest`  

---

### **Installation**

1. **Clone the repository**:
   ```bash
   git clone https://github.com/caterer-z-t/SOM_Seq_Sim.git
   cd /path/to/SOM_Seq_Sim
   ```

2. **Create a virtual environment (optional but recommended):**

    If you're using `conda`, `mamba`, or `micromamba`, you can create the environment directly from the `env.yml` file:
    ``` bash 
    conda env create -f env.yml
    ```

3. Activate environment using conda (please change {conda} to be whatever used to create the environemnt):
    ``` bash
    conda activate som-sim-env
    ```

### **Generating Sequencing Data**
  > Note: The code (`Seq_Sim/seq_sim.py` and `Seq_Sim/utils/seq_sim_utils.py`) is based on R scripts from the [Zhang Lab](https://fanzhanglab.org/). These R scripts have been modified and streamlined for this project. The biological relevance may not be fully retained, and it serves as a showcase for potential sequencing simulations. For more information please see [`Seq_Sim/README.md`](Seq_Sim/README.md)
  
You can generate simulated sequencing data by running the following command:

``` bash
/Seq_Sim/simulations.sh -c /Seq_Sim/config.yml
```
Alternatively, you can run:
``` bash
python /Seq_Sim/seq_sim.py \ 
      --num_samples 30                   \ # num samples 
      --fold_change 0.5                  \ # fold change between disease and healthy samples 
      --config_file /Seq_Sim/config.yml    # configuration file
```

By default, the output CSV files will be saved in the `SOM_Seq_Sim/data/` directory.

### Running SOMs
For running Self-Organizing Maps, whether on the generated sequencing data or your own data, please consult the SOM directory's [README file](SOM/).

## Documentation
For our documentation, please see the [readthedocs](https://som-seq-sim.readthedocs.io/en/latest/) page.

## Contributing 
To contribute to this project, please see the [CONTRIBUTING.md](CONTRIBUTING.md) file.

## Code of Conduct
For our code of conduct for this project, please see the [CODE_OF_CONDUCt.md](CODE_OF_CONDUCT.md) file. 

## Testing
Run the following command to execute tests:
```bash
pytest test/ -v
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
      <b>Con</b>
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

## Project Structure

``` bash
SOM_Seq_Sim/
├── docs/
│   ├── img/                              # Contains images used in Read the Docs
│   ├── pipelines/                        # Additional MD files for pipelines and packages
│   ├── index.md                          # Initial page for Read the Docs
│   └── requirements.txt                  # Dependencies for Read the Docs 
├── SOM/
│   ├── utils/
│   │       └── som_utils.py              # Contains the core implementation of the `SOM` class
│   ├── data/              
│   │   ├── iris_categorical_data.csv     # Categorical data for the Iris dataset
│   │   ├── iris_training_data.csv        # Training data for the Iris dataset
│   │   ├── seq_sim_categorical_data.csv  # Categorical data for the sequencing dataset
│   │   ├── seq_sim_training_data.csv     # Training data for the sequencing dataset
│   │   ├── titanic_categorical_data.csv  # Categorical data for the Titanic dataset
│   │   └── titanic_training_data.csv     # Training data for the Titanic dataset
│   ├── examples/           
│   │   ├── output/
│   │   │   ├── iris/                     # Directory containing output figures from `som_example_iris.py`
│   │   │   ├── seq/                      # Directory containing output figures from `som_example_seq.py`
│   │   │   └── titanic/                  # Directory containing output figures from `som_example_titanic.py`
│   │   ├── som_example_iris.py           # Example script for applying the SOM to the Iris dataset
│   │   ├── som_example_seq.py            # Example script for applying the SOM to the Sequencing dataset
│   │   └── som_example_titanic.py        # Example script for applying the SOM to the Titanic dataset
│   ├── README.md/                        # README file for SOM directory
│   └── som.py                            # A command-line interface for fitting and visualizing SOMs. 
├── Seq_Sim/
│   ├── utils/                            # folder containg additional utility functions
│   │   └── seq_sim_utils.py              # contains all necessary functions for simulation.R
│   ├── config.yml                        # configuration file use to specify all components in simulation.R
│   ├── README.md                         # information specific to running this simulation sequencing generation code
│   ├── seq_sim.py                        # main function which runs the seqeucning generation code
│   └── simulation.sh                     # shell script for running sequencing generating code multiple times, see Seq_Sim/README.md for more information
├── test/
│   ├── Seq_Sim/ 
│   │   ├── test_seq_sim_utils.py
│   │   └── test_seq_sim.py
│   └── SOM/
│       └── test_SOM.py                   # Tests for SOM functionality
├── .gitignore                            # Git ignore file to exclude files and directories from Git
├── .readthedocs.yaml                     # YAML file specifying dependencies for Read the Docs
├── CODE_OF_CONDUCT.md                    # markdown file explaining our code of conduct
├── CONTRIBUTING.md                       # markdown file explaining how to contribute to this repo
├── env.yml                               # YAML file specifying dependencies for the entire package
├── README.md                             # README file
└── LICENSE                               # License for the project
```

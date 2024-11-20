<img src="assets/som_package_logo.png" alt="SOM SIM SEQ Logo" width="25%">

# **Self-Organizing Maps for Genetic Sequencing Simulation (SOM-Seq)**

**Simulation of Genetic Sequencing Data using Statistically Similar Distributions from Original Datasets**

[![GitHub release](https://img.shields.io/github/v/release/caterer-z-t/SOM_Seq_Sim)](https://github.com/caterer-z-t/SOM_Seq_Sim/releases) 
[![GitHub Release Date](https://img.shields.io/github/release-date/caterer-z-t/SOM_Seq_Sim)](https://github.com/caterer-z-t/SOM_Seq_Sim/releases) 
[![Tests](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/test.yml?label=tests&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/actions)

---

## **Course:**
[CSCI 6118](https://github.com/swe4s): Software Engineering for Scientists at the [University of Colorado Boulder](https://www.colorado.edu/)

---

## **Project Overview**

Welcome to **SOM-Seq**, a cutting-edge approach to **simulating genetic sequencing data** using **Self-Organizing Maps (SOMs)**! 🎉  
This project leverages SOMs to capture the statistical distributions of genetic data and generate **new, statistically similar datasets**. Whether you're diving into **population genetics**, **personalized medicine**, or **genomic data analysis**, SOM-Seq can help create synthetic, yet realistic datasets, a critical tool when real-world data is scarce. 🚀

---

## **Table of Contents**
- [Project Setup](#project-setup)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Generating Sequencing Data](#generating-sequencing-data)
- [Running SOMs](#running-soms)
- [Contributing](#contributing)
- [Code of Conduct](#code-of-conduct)
- [Contributors](#contributors)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## **Usage**

### **Project Setup**

#### **Requirements**

This project depends on several libraries and software packages. You can find the list in the `env.yml` file, which makes setting up the environment a breeze.

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
### **Usage**

#### Generating Sequencing Data
  > Note: The code (`Seq_Sim/seq_sim.py` and `Seq_Sim/utils/seq_sim_utils.py`) is based on R scripts from the [Zhang Lab](https://fanzhanglab.org/). These R scripts have been modified and streamlined for this project. The biological relevance may not be fully retained, and it serves as a showcase for potential sequencing simulations. For more information please see [`Seq_Sim/README.md`](Seq_Sim/README.md)

  **Step 1. Generate Simulated Single Cell Sequencing Data**
  
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

#### Running SOMs
For running Self-Organizing Maps on the generated sequencing data, further instructions are required. Please refer to the documentation or instructions provided by @madrpernat

## Contributing 
To contribute to this project, please see the [CONTRIBUTING.md](CONTRIBUTING.md) file.

## Code of Conduct
For our code of conduct for this project, please see the [CODE_OF_CONDUCt.md](CODE_OF_CONDUCT.md) file. 


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


### Project Structure
``` bash 
SOM_Seq_Sim/
├── docs/
│   ├── img/                     # Contains images used in Read the Docs
│   ├── pipelines/               # Additional MD files for pipelines and packages
│   ├── index.md                 # Initial page for Read the Docs
│   └── requirements.txt         # Dependencies for Read the Docs 
├── SOM/
├── Seq_Sim/
│   ├── utils/                  # folder containg additional utility functions
│   │   └── seq_sim_utils.py    # contains all necessary functions for simulation.R
│   ├── config.yml              # configuration file use to specify all components in simulation.R
│   ├── README.md               # information specific to running this simulation sequencing generation code
│   ├── seq_sim.py              # main function which runs the seqeucning generation code
│   └── simulation.sh           # shell script for running sequencing generating code multiple times, see Seq_Sim/README.md for more information
├── test/
│   ├── Seq_Sim/ 
│   │   ├── test_seq_sim_utils.py
│   │   └── test_seq_sim.py 
├── .gitignore                  # Git ignore file to exclude files and directories from Git
├── .readthedocs.yaml           # YAML file specifying dependencies for Read the Docs
├── CODE_OF_CONDUCT.md          # markdown file explaining our code of conduct
├── CONTRIBUTING.md             # markdown file explaining how to contribute to this repo
├── env.yml                     # YAML file specifying dependencies for the entire package
├── README.md                   # README file
└── LICENSE                     # License for the project
```
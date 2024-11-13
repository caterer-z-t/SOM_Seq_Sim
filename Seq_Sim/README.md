# Sequence Simulation

here lies the directory containing the simulation sequencing code. This code is designed to generate sequencing data which can be implemented in various forms of stasticial and mathematical purposes such as in the development of novel software for analysis of single cell RNA sequencing datasets in which you want the ideal situation to develop and employ your code. additionaly, in the rhelm of comparing many forms of the same type of analysis. here having predefined relationships in the inhernet dataset where the goal is to find them not necessarily who is able to find the move obsucre relationship. Subsequently, this code is designed to depoly that exact concept. Below, we describe the process to generate your own sequencing data.

1. ***Install python, and R***

please reference this section to install python and R respectively. 

2. ***Install necessary R packages***

Since this section of our code is written in R, please install the necessary R packages through the conda environment which is required to run the rest of this software. These packages should have been installed directly using the `env.yml`. To double check these packages are installed and if not, please install the following 

``` bash
cd /path/to/SOM_Seq_Sim/Seq_Sim/
R
```

``` R 
> library(reticulate) #required for the conda environment
> # if not installed then run
> # install.packages('reticulate')

> use_condaenv( " # name of your working environment ",
                required = TRUE)

> # Load the libraries and install any packages if needed
> suppressWarnings({
    library(dplyr)
    library(tidyr)
    library(tidyverse)
    library(ggalluvial)
    library(ggrepel)
    library(MASS)
    library(caret)
    library(Seurat)
    library(ggplot2)
    library(glue)
    library(stevemisc)
    library(stevedata)
    library(lme4)
    library(broom.mixed)
    library(moments)
    library(ggpubr)
    library(ggh4x)
    library(doParallel)
    library(pbapply)
    library(variancePartition)
    library(pheatmap)
    library(purrr)
    library(data.table)
    library(presto)
    library(harmony)
    library(foreach)
})

> quit() # to exit the R interactive environment in the terminal
```
:tip: if it is needed please complete this using R studio for easier installational issue naviagation. 

3. ***Modify the `config.yml`***

Once the working environment is set up and ready with all necessary (R, specifically) packages installed. We will want to modify the `config.yml` file. Please see this file to understand its structure and contents. Please copy this file, and modify it accordingly to your specific use case of this software. It is important to note that there are 2 parameters of importance

``` yaml
num_samples:
  - 10
  - 20
  - 30
fold_changes:
  - 0.1
  - 0.75
  - 1.5
  - 3
```

`num_samples` and `fold_changes`, these parameters will specify how many files in total will be generated. In this case we will have 12 seperate experiments (`fold_change` * `num_samples`) and each experiement produces 4 files (this can be modified if desired, will explain later) so in total there will be 48 files generated in total. Please change this accordingly. 

Once you are satisified with the yaml file, please save it accordingly and name it appropiately. 

4. ***Running of `simulation.R` and `simulation.sh`***

Now we have a properly configured yaml file we can run this script to generate our sequencing data in 2 ways. 

A. directly, through bash

B. through the `simulation.sh` shell script. 

    
  A. For running `simulatiom.R` directly through bash, please use the following. Please note that this method will run the script for one time resulting in 4 saved csv files, regardless of how many `fold_changes` or `num_samples` are in the `config.yaml` file, this is because these are arguments for running the R script. 

    ``` bash
    # requirements for running this script
    # r_script_file, file path to the R script, in this case it is the file path to `simulation.R`
    # num_samp, this is the number of samples per patient in the analysis
    # fold_change, this is the fold change between the healthy patients and the diseased patients
    # config_file, this is the file path to the `config.yml` file

    Rscript /filepath/to/simulation.R num_samp fold_change /filepath/to/config.yml
    ```
  B. For running `simulation.sh` please use the following, this will require 2 argument. 1. the filepath to the `simulation.sh` file, and 2. the filepath to the `config.yml` file. Please use the following command to run this script. 

  ``` bash
  /filepath/to/simulation.sh /filepath/to/config.yml
  ```



Rscript /Users/zc/Library/CloudStorage/OneDrive-UCB-O365/academics/2024_fall/CSCI_6118/SOM_Seq_Sim/Seq_Sim/simulation.R 10 3 /Users/zc/Library/CloudStorage/OneDrive-UCB-O365/academics/2024_fall/CSCI_6118/SOM_Seq_Sim/Seq_Sim/config.yml
# Self-Organizing Maps (SOM) Framework

This repository provides an implementation of Self-Organizing Maps (SOMs), a clustering and visualization tool for high-dimensional data. The framework includes utilities, example scripts, and sample data to demonstrate the SOM's capabilities.

---

## Directory Structure

``` bash
SOM/
├── utils/
│   └── som_utils.py 
├── data/              
│   ├── iris_categorical_data.csv          
│   ├── iris_training_data.csv     
│   ├── seq_sim_categorical_data.csv
│   ├── seq_sim_training_data.csv 
│   ├── titanic_categorical_data.csv
│   └── titanic_training_data.csv 
├── examples/           
│   ├── output/
│   │   ├── iris/
│   │   ├── seq/
│   │   └── titanic/
│   ├── som_example_iris.py
│   ├── som_example_seq.py
│   └── som_example_titanic.py             
└── som.py                  
```

## Features

- **Core Implementation**: The `som_utils.py` file contains the `SOM` class, which provides functionalities for training, scaling, and visualizing Self-Organizing Maps.
- **Command-Line Interface**: The `som.py` script enables running SOM operations from the command line for integration into pipelines or automation.
- **Example Scripts**: The `examples` directory provides ready-to-run scripts to demonstrate training, clustering, and visualization using the SOM framework.
- **Data Support**: The `data` directory includes sample datasets to test and validate the functionality of the SOM framework.

---
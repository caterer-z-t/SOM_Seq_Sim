Creating SOMs
=============

This code can be utilized in two main ways, depending on your goals. The first is through a command-line interface (CLI), which offers a streamlined and efficient approach for training and visualizing a SOM, as well as performing hyperparameter tuning. The second option involves writing a custom script that leverages the ``SOM`` class in ``som_utils.py``. This method is more flexible and is recommended if you plan to perform additional analyses beyond hyperparameter tuning, clustering, and visualization—for example, implementing secondary clustering of the SOM, which is outside the scope of this project.

Using the Command-Line Interface
--------------------------------

Assuming you are located in the ``SOM_Seq_Sim`` repository, navigate to the SOM directory:

.. code-block:: bash

   cd SOM

Within this directory, the ``som.py`` script provides a command-line interface (CLI) for training and visualizing SOMs. The CLI supports both direct SOM fitting with user-specified hyperparameters **or** tuning to identify the optimal hyperparameter configuration for your dataset.

### Inputs

#### Required Inputs
- **Training Data** (``-t`` or ``--train_dat``):  
  Path to a CSV file containing numerical training data. Each row represents a data point, and each column represents a feature.

- **Output Directory** (``-o`` or ``--output_directory``):  
  Path to save SOM results, including metrics and visualizations.

- **Scaling Method** (``-s`` or ``--scale_method``):  
  Scaling method for the data. Accepted values: ``zscore`` or ``minmax``.

- **Grid Dimensions** (``-x`` and ``-y`` or ``--x_dim`` and ``--y_dim``):  
  Number of neurons in the x and y dimensions of the SOM grid.

- **Topology** (``-p`` or ``--topology``):  
  Grid topology for the SOM. Accepted values: ``rectangular`` or ``hexagonal``.

- **Neighborhood Function** (``-n`` or ``--neighborhood_fnc``):  
  Neighborhood function for training the SOM. Accepted values: ``gaussian`` or ``bubble``.

- **Epochs** (``-e`` or ``--epochs``):  
  Number of training epochs.

#### Optional Inputs
- **Categorical Data** (``-c`` or ``--other_dat``):  
  Path to a CSV file containing categorical metadata associated with the training data. These categories are not used in training but can be visualized on the SOM grid. This file must have the same number of rows as in Training Data.

- **Component Plane Plots** (``-m`` or ``--plot_component_planes``):  
  Whether to generate component plane plots for each feature. Defaults to ``True``.

### Usage

#### Hyperparameter Tuning
If multiple values are provided for any hyperparameter, the CLI will perform hyperparameter tuning by testing all combinations of the specified values.

.. code-block:: bash

   python som.py -t data/iris_training_data.csv \
                 -c data/iris_categorical_data.csv \
                 -o examples/output/iris/ \
                 -s zscore minmax \
                 -x 3 5 7 \
                 -y 2 4 \
                 -p rectangular hexagonal \
                 -n gaussian bubble \
                 -e 50 100

This command evaluates all combinations of the specified hyperparameters, identifies the best configuration based on PVE and Topographic Error, and retrains the SOM using those parameters.

#### Direct SOM Fitting
If a single value is provided for each hyperparameter, the CLI will directly fit the SOM with those parameters.

.. code-block:: bash

   python som.py -t data/iris_training_data.csv \
                 -c data/iris_categorical_data.csv \
                 -o examples/output/iris/ \
                 -s zscore \
                 -x 7 \
                 -y 4 \
                 -p hexagonal \
                 -n gaussian \
                 -e 50

### Outputs

#### Metrics
- **Percent Variance Explained (PVE)**: Indicates how well the SOM represents the variance in the dataset. This is comparable to Within Cluster Sum of Squares used in other clustering methods.
- **Topographic Error**: Measures how well the SOM preserves the topology of the input data. This is calculated as the percentage of data points whose second closest neuron is not adjacent to the neuron it's mapped to.
- **Note**: If hyperparameter tuning is performed, these metrics will be output for each configuration tested.

#### Visualizations
- **Component Plane Plots**: Figures displaying the distribution of each feature across the SOM grid as heatmaps.
- **Categorical Data Overlays**: Visualizations showing the distribution of external categories (if provided) across the SOM grid.
- **Note**: If hyperparameter tuning is performed, these visualizations will be generated only for the best hyperparameter configuration.

Writing Your Own Script
-----------------------

If you prefer greater flexibility in how you utilize the SOM, you can write your own script that imports the ``SOM`` class from ``som_utils.py``. This approach allows you to fully customize the SOM workflow and incorporate additional analyses or operations beyond the scope of the CLI. For instance, you could perform secondary clustering, experiment with alternative visualization methods, or integrate the SOM with other tools in your data analysis pipeline.

Refer to the :doc:`examples <examples/>` directory for sample scripts, which demonstrate how to set up and use the ``SOM`` class for various datasets.

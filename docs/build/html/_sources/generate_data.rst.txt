Generating Sequencing Data
==========================

.. note::

   The code (``Seq_Sim/seq_sim.py`` and ``Seq_Sim/utils/seq_sim_utils.py``) is based on R scripts from the `Zhang Lab <https://fanzhanglab.org/>`_. These R scripts have been modified and streamlined for this project. The biological relevance may not be fully retained, and it serves as a showcase for potential sequencing simulations. For more information, please see :doc:`Seq_Sim/README.md <Seq_Sim/README.md>`.

You can generate simulated sequencing data by running the following command:

.. code-block:: bash

   /Seq_Sim/simulations.sh -c /Seq_Sim/config.yml

Alternatively, you can run:

.. code-block:: bash

   python /Seq_Sim/seq_sim.py \ 
         --num_samples 30                   \ # num samples 
         --fold_change 0.5                  \ # fold change between disease and healthy samples 
         --config_file /Seq_Sim/config.yml    # configuration file

By default, the output CSV files will be saved in the ``SOM_Seq_Sim/data/`` directory.

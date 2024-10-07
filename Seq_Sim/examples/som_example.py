import pandas as pd
from Seq_Sim.utils.SOM import SOM

train_dat = pd.read_csv("./data/simulation/simulated_pc_NAM_HARMONIZED_jade_woi_advanced_numCell100_numSamp10_sd0.1_fci0.1_clusterRatio0.7_diseaseRatio0.05.csv")
other_dat = pd.read_csv("./data/simulation/simulated_umap_data_jade_woi_advanced_numCell100_numSamp10_sd0.1_fci0.1_clusterRatio0.7_diseaseRatio0.05.csv")

som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method="zscore",
    x_dim=8,
    y_dim=6,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=12
)

som.train_map()



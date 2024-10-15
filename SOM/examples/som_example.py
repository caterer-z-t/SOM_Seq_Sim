import pandas as pd
from SOM.utils.SOM import SOM

# Training data: Only considering two features here. This data has already undergone PCA and so 
# there aren't really any patterns between the features. This is more just a proof of concept
# that the training and plotting works.
train_dat = pd.read_csv("./data/simulation/simulated_pc_NAM_HARMONIZED_jade_woi_advanced_numCell100_numSamp10_sd0.1_fci0.1_clusterRatio0.7_diseaseRatio0.05.csv")
train_dat = train_dat.iloc[:, 0:2] 

# Other data: Data that is NOT used for training but still could be interesting to project onto the
# map. No plotting methods for this yet.
other_dat = pd.read_csv("./data/simulation/simulated_umap_data_jade_woi_advanced_numCell100_numSamp10_sd0.1_fci0.1_clusterRatio0.7_diseaseRatio0.05.csv")

som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method="zscore",
    x_dim=8,
    y_dim=6,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=5
)

som.train_map()

som.plot_component_planes(
    output_dir="Seq_Sim/examples/output_figs"
)

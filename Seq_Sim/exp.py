from Seq_Sim.utils.TimeSeriesSimSeq import TimeSeriesSimSeq
from Seq_Sim.utils.SimSeq import SimSeq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Example usage
data = {
    "numerical_col": np.random.randn(10),
    "categorical_col": ["A", "B", "A", "C", "B", "A", "C", "C", "B", "A"],
    "group": [
        "healthy",
        "disease",
        "healthy",
        "healthy",
        "disease",
        "healthy",
        "disease",
        "disease",
        "healthy",
        "disease",
    ],
}
df = pd.DataFrame(data)

# Initialize the TimeSeriesSimSeq object
sim_seq = TimeSeriesSimSeq(df)

# Generate synthetic data
num_rows = 5000
num_time_points = 5
synthetic_data_list = sim_seq.generate_synthetic_time_series(
    num_rows=num_rows, num_time_points=num_time_points, group_col="group"
)


# Initialize the SimSeq object for plotting
# Assuming you want to plot comparisons of the last time point in synthetic_data_list
# and compare it with the original DataFrame
sim_seq_plotter = SimSeq(df, synthetic_df=synthetic_data_list[-1])

# Plot the original and synthetic data
sim_seq_plotter.plot_numerical_distributors()
sim_seq_plotter.plot_categorical_distributions()
sim_seq_plotter.plot_comparisons()

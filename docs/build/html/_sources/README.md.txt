<img src="assets/som_package_logo.png" alt="SOM-Seq Logo" width="20%">

# SOM-Seq

[![GitHub release](https://img.shields.io/github/v/release/caterer-z-t/SOM_Seq_Sim)](https://github.com/caterer-z-t/SOM_Seq_Sim/releases)
[![Tests](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/test.yml?label=tests&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/actions)
[![Build](https://img.shields.io/github/actions/workflow/status/caterer-z-t/SOM_Seq_Sim/build.yml?label=build&logo=github)](https://github.com/caterer-z-t/SOM_Seq_Sim/actions)
[![codecov](https://codecov.io/gh/caterer-z-t/SOM_Seq_Sim/branch/main/graph/badge.svg)](https://codecov.io/gh/caterer-z-t/SOM_Seq_Sim)
[![License](https://img.shields.io/github/license/caterer-z-t/SOM_Seq_Sim)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-online-blue?logo=readthedocs)](https://som-seq-sim.readthedocs.io)

**SOM-Seq** is a Python toolbox that combines synthetic single-cell sequencing data generation (`Seq_Sim`) with Self-Organizing Map (SOM) clustering and visualization (`SOM`). The two modules can be used together or independently.

---

## Installation

**From source:**
```bash
git clone https://github.com/caterer-z-t/SOM_Seq_Sim.git
cd SOM_Seq_Sim
pip install .
```

---

## Quick Start

### Generate Sequencing Data

```bash
seq-sim --num_samples 30 --fold_change 0.5 --config_file Seq_Sim/config.yml
```

Or equivalently:
```bash
python Seq_Sim/seq_sim.py \
    --num_samples 30 \
    --fold_change 0.5 \
    --config_file Seq_Sim/config.yml
```

Output CSV files are written to the directory specified in `config.yml` (default: `data/`).

See [`Seq_Sim/config.yml`](Seq_Sim/config.yml) for all configurable parameters (cell-type counts, batch structure, disease proportions, number of features, etc.).

---

### Fit a SOM

**Single hyperparameter set:**
```bash
som-fit \
    -t data/seq_sim_training_data.csv \
    -c data/seq_sim_categorical_data.csv \
    -o output/ \
    -s zscore -x 5 -y 4 -p hexagonal -n gaussian -e 100
```

**Hyperparameter tuning (pass multiple values; best combination is selected automatically):**
```bash
som-fit \
    -t data/seq_sim_training_data.csv \
    -c data/seq_sim_categorical_data.csv \
    -o output/ \
    -s zscore minmax -x 3 5 7 -y 3 5 -p rectangular hexagonal -n gaussian -e 50 100
```

Or equivalently use `python SOM/som.py` with the same flags.

**Python API:**
```python
import pandas as pd
from SOM.utils.som_utils import SOM

train = pd.read_csv("data/seq_sim_training_data.csv")
meta  = pd.read_csv("data/seq_sim_categorical_data.csv")

som = SOM(
    train_dat=train,
    other_dat=meta,
    scale_method="zscore",
    x_dim=5,
    y_dim=4,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=100,
)
som.train_map()

print(f"PVE:               {som.calculate_percent_variance_explained():.1f}%")
print(f"Topographic error: {som.calculate_topographic_error():.3f}")

som.plot_component_planes(output_dir="output/")
som.plot_categorical_data(output_dir="output/")
```

---

## CLI Reference

### `seq-sim` / `python Seq_Sim/seq_sim.py`

| Flag | Description |
|------|-------------|
| `--num_samples` | Number of subjects to simulate |
| `--fold_change` | Disease-associated fold change magnitude |
| `--config_file` | Path to `config.yml` |

### `som-fit` / `python SOM/som.py`

| Flag | Description |
|------|-------------|
| `-t` | Path to training data CSV (numeric features) |
| `-c` | Path to categorical/metadata CSV (optional) |
| `-o` | Output directory for plots and metrics |
| `-s` | Scaling method: `zscore` or `minmax` (one or more) |
| `-x` | SOM x-dimension (one or more integers) |
| `-y` | SOM y-dimension (one or more integers) |
| `-p` | Topology: `rectangular` or `hexagonal` (one or more) |
| `-n` | Neighborhood function: `gaussian` or `bubble` (one or more) |
| `-e` | Number of training epochs (one or more integers) |
| `-m` | Generate component plane plots (default: `True`) |

When more than one value is provided for any hyperparameter, the CLI performs a grid search and selects the best combination by `PVE − 100 × topographic_error`.

---

## Documentation

Full API reference and tutorials: [som-seq-sim.readthedocs.io](https://som-seq-sim.readthedocs.io)

---

## Testing

```bash
pytest test/ -v
```

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please follow the [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Citation

If you use SOM-Seq in your research, please cite:

> Caterer Z., Pernat M., Hurd V. (2024). *SOM-Seq: A Python Toolbox for Single-Cell Sequencing Simulation and Self-Organizing Map Analysis*. <!-- TODO: add journal/DOI after JOSS acceptance -->

A machine-readable citation is available in [CITATION.cff](CITATION.cff).

---

## License

MIT — see [LICENSE](LICENSE).

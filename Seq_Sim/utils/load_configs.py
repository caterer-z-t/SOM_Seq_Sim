###############################################################
###
###                     Imports
###
###############################################################


import argparse
import yaml


###############################################################
###
###                     Function
###
###############################################################


def load_config(args):
    """Load configuration either from YAML or command-line arguments."""
    if args.config:
        # Load the YAML configuration file
        with open(args.config, "r") as file:
            config = yaml.safe_load(file)

        # Extract configuration options from the YAML file
        input_path = config.get("input")
        output_path = config.get("output")
        rows = config.get("rows")
        output_dir = config.get("output_dir")
    else:
        # Use command-line arguments as configuration
        input_path = args.input
        output_path = args.output
        rows = args.rows
        output_dir = args.output_dir

    return input_path, output_path, rows, output_dir

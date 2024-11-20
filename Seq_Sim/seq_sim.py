import sys

from Seq_Sim.utils.seq_sim_utils import load_config, validate_arguments, generate_and_save_features, arg_parser

def main() -> None:
    """Main function to generate and save features for the given number of samples and fold change.
    """
    args = arg_parser()
    try:
        num_samples, fold_change, config_file = validate_arguments(list(args.config, args.num_samples, args.fold_change))
        config = load_config(config_file)
        generate_and_save_features(num_samples, fold_change, config)
        print("Data generation and saving completed successfully.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

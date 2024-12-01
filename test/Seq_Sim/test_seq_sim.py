import unittest
from unittest.mock import patch, MagicMock
from argparse import Namespace
import sys

# Assuming your script is in a file called `seq_sim_main.py`:
from Seq_Sim.seq_sim import main


class TestSeqSimScript(unittest.TestCase):

    @patch("Seq_Sim.utils.seq_sim_utils.arg_parser")
    @patch("Seq_Sim.utils.seq_sim_utils.generate_and_save_features")
    @patch("Seq_Sim.utils.seq_sim_utils.load_config")
    @patch("Seq_Sim.utils.seq_sim_utils.validate_arguments")
    def test_main_all_paths(
        self,
        mock_validate_arguments,
        mock_load_config,
        mock_generate_and_save_features,
        mock_arg_parser,
    ):
        # Mock the parsed arguments
        mock_arg_parser.return_value = Namespace(
            config_file="Seq_Sim/config.yaml", num_samples=10, fold_change=2.0
        )

        # Mock return values for other functions
        mock_validate_arguments.return_value = (10, 2.0, "Seq_Sim/config.yaml")
        mock_load_config.return_value = {"dummy_key": "dummy_value"}
        mock_generate_and_save_features.return_value = None

        # Simulate command-line arguments
        sys.argv = [
            "Seq_Sim/seq_sim.py",
            "--config_file",
            "Seq_Sim/config.yaml",
            "--num_samples",
            "10",
            "--fold_change",
            "2.0",
        ]

    @patch("Seq_Sim.utils.seq_sim_utils.arg_parser")
    @patch("Seq_Sim.utils.seq_sim_utils.validate_arguments")
    def test_main_with_exception(self, mock_validate_arguments, mock_arg_parser):
        # Simulate argument parser output
        mock_arg_parser.return_value = MagicMock(
            config="config.yaml", num_samples=10, fold_change=2.0
        )

        # Make validate_arguments raise an exception
        mock_validate_arguments.side_effect = Exception("Invalid arguments")

        # Run main and check exception handling
        sys.argv = [
            "Seq_Sim/seq_sim.py",
            "--config",
            "Seq_Sim/config.yaml",
            "--num_samples",
            "10",
            "--fold_change",
            "2.0",
        ]
        with self.assertRaises(SystemExit) as cm:
            main()
        self.assertEqual(cm.exception.code, 1)

    @patch("Seq_Sim.utils.seq_sim_utils.arg_parser")
    @patch("Seq_Sim.utils.seq_sim_utils.load_config")
    def test_main_load_config_failure(self, mock_load_config, mock_arg_parser):
        # Simulate argument parser output
        mock_arg_parser.return_value = MagicMock(
            config="config.yaml", num_samples=10, fold_change=2.0
        )

        # Make load_config raise an exception
        mock_load_config.side_effect = Exception("Failed to load config")

        # Run main and check exception handling
        sys.argv = [
            "Seq_Sim/seq_sim.py",
            "--config",
            "Seq_Sim/config.yaml",
            "--num_samples",
            "10",
            "--fold_change",
            "2.0",
        ]
        with self.assertRaises(SystemExit) as cm:
            main()
        self.assertEqual(cm.exception.code, 1)


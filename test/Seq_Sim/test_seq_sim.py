import unittest
from unittest.mock import patch
import sys

# Assuming your script is in a file called `seq_sim_main.py`:
from Seq_Sim.seq_sim import main

class TestSeqSimScript(unittest.TestCase):

    @patch('Seq_Sim.utils.seq_sim_utils.load_config')
    @patch('Seq_Sim.utils.seq_sim_utils.validate_arguments')
    def test_main_with_exception(self, mock_validate_arguments, _):
        # Simulating a failure in the argument validation
        sys.argv = ['seq_sim_main.py', 'config.yaml']

        # Make `validate_arguments` raise an exception
        mock_validate_arguments.side_effect = Exception("Invalid arguments")

        # Call the main function and test if exception handling works
        with self.assertRaises(SystemExit):
            main()

if __name__ == '__main__':
    unittest.main()

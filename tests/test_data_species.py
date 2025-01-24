import unittest
from unittest.mock import MagicMock, patch
import numpy as np
from melektrodica.collector import DataSpecies


class TestDataSpecies(unittest.TestCase):
    """
    Unit tests for the DataSpecies class to verify initialization, data processing,
    and error handling based on various scenarios.
    """

    def setUp(self):
        """
        Initialize resources for the test cases:
        - Mock parameters and writer objects.
        - Simulate a species file path to be used in the tests.
        """
        self.mock_writer = MagicMock()  # Mock for logging or writing messages
        self.mock_parameters = MagicMock()  # Mock for input parameters
        self.mock_parameters.thermochemical = True  # Enable thermochemical calculations by default
        self.mock_parameters.g_formation = True  # Assume energy of formation is provided
        self.mock_species_file = "test_species_file.md"  # Simulated input file for testing

    @patch("melektrodica.collector.Collector.raw_data")
    def test_initialization_valid_md_file(self, mock_raw_data):
        """
        Test the valid initialization of the DataSpecies object
        using a properly structured Markdown file with species data.
        """
        # Simulate the data returned by the Collector.raw_data function
        mock_raw_data.return_value = (
            ["Species", "RPACe", "DG_formation", "c0", "Catalyst"],
            np.array([
                ["O2", "R", 0.0, 0.5307, ""],
                ["H+", "R", 0.0, 1.0, ""],
                ["H2O", "P", 0.0, 1.0, ""],
                ["O*", "A", -0.343, "", "Pt"],
                ["OH*", "A", -0.376, "", "Pt"],
                ["Pt", "C", "", "", ""],
                ["e-", "e", "", "", ""]
            ])
        )

        # Initialize the DataSpecies object
        data_species = DataSpecies(self.mock_species_file, self.mock_parameters, self.mock_writer)

        # Assert the categorized species lists
        self.assertEqual(data_species.reactants, ["O2", "H+"])
        self.assertEqual(data_species.products, ["H2O"])
        self.assertEqual(data_species.adsorbed, ["O*", "OH*"])
        self.assertEqual(data_species.catalyst, ["Pt"])

        # Assert initial concentrations for reactants and products
        np.testing.assert_array_equal(data_species.c0_reactants, np.array([0.5307, 1.0]))
        np.testing.assert_array_equal(data_species.c0_products, np.array([1.0]))

        # Assert the full list of species in the system
        self.assertEqual(data_species.list, ["O2", "H+", "H2O", "O*", "OH*", "Pt", "e-"])

        # Assert Gibbs free energy values for reactants, adsorbed species, and products
        np.testing.assert_array_equal(data_species.g_formation_rct, np.array([0.0, 0.0]))
        np.testing.assert_array_equal(data_species.g_formation_ads, np.array([-0.343, -0.376]))
        np.testing.assert_array_equal(data_species.g_formation_prd, np.array([0.0]))

        # Assert logging calls to the writer mock
        self.mock_writer.message.assert_any_call(f"Reading species data from file: {self.mock_species_file}")

    @patch("melektrodica.collector.Collector.raw_data")
    def test_initialization_missing_columns(self, mock_raw_data):
        """
        Test initialization when the input file is missing required columns.
        This should raise a ValueError.
        """
        # Simulate missing required columns in the input data
        mock_raw_data.return_value = (
            ["Species", "RPACe"],
            np.array([
                ["O2", "R"],
                ["H+", "R"]
            ])
        )

        # Expect a ValueError due to missing columns
        with self.assertRaises(ValueError) as context:
            DataSpecies(self.mock_species_file, self.mock_parameters, self.mock_writer)

        self.assertIn(
            "Missing required column 'c0' in the header of file 'test_species_file.md'.",
            str(context.exception)
        )

    @patch("melektrodica.collector.Collector.raw_data")
    def test_initialization_no_thermochemical(self, mock_raw_data):
        """
        Test initialization when thermochemical calculations are disabled
        (`thermochemical = False`).
        """
        # Disable thermochemical processing
        self.mock_parameters.thermochemical = False

        # Simulate valid input data without thermochemical parameters
        mock_raw_data.return_value = (
            ["Species", "RPACe", "c0", "Catalyst"],
            np.array([
                ["O2", "R", 0.5307, ""],
                ["H+", "R", 1.0, ""],
                ["H2O", "P", 1.0, ""],
            ])
        )

        # Initialize the DataSpecies object
        data_species = DataSpecies(self.mock_species_file, self.mock_parameters, self.mock_writer)

        # Assert that thermochemical data is not used
        self.assertEqual(data_species.reactants, ["O2", "H+"])
        self.assertEqual(data_species.products, ["H2O"])
        self.assertEqual(data_species.adsorbed, [])
        self.assertEqual(data_species.catalyst, [])
        np.testing.assert_array_equal(data_species.c0_reactants, np.array([0.5307, 1.0]))
        np.testing.assert_array_equal(data_species.c0_products, np.array([1.0]))
        self.assertEqual(data_species.list, ["O2", "H+", "H2O", "e-"])

    @patch("melektrodica.collector.Collector.raw_data")
    def test_initialization_invalid_file(self, mock_raw_data):
        """
        Test initialization when the input species file is invalid or missing.
        This should raise a FileNotFoundError.
        """
        # Simulate that the file is not found
        mock_raw_data.side_effect = FileNotFoundError("Species file not found")

        # Expect a FileNotFoundError during initialization
        with self.assertRaises(FileNotFoundError) as context:
            DataSpecies(self.mock_species_file, self.mock_parameters, self.mock_writer)

        self.assertIn("Species file not found", str(context.exception))

if __name__ == "__main__":
    unittest.main()

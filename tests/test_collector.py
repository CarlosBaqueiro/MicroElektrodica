import numpy as np
import os
import unittest
from unittest.mock import MagicMock, patch
from melektrodica import *


class TestCollector(unittest.TestCase):
    """
    Unit test class for Collector to validate the functionality of its methods
    including normal use cases, edge cases, and exception scenarios.
    """

    def setUp(self):
        self.test_directory = "test_data"
        self.writer = MagicMock()

        # Mock for header and raw data
        self.header = ["Column1", "Column2", "Column3"]
        self.raw_data = np.array([[1, 2, 3], [4, 5, 6]])

        # Mock file content
        self.sample_file_content = [
            "| Column1 | Column2 | Column3 |\n",
            "|---|---|---|\n",
            "| 1 | 2 | 3 |\n",
            "| 4 | 5 | 6 |\n",
        ]

        # Create test directory if not exists
        if not os.path.exists(self.test_directory):
            os.makedirs(self.test_directory)

        # Create fake files required by Collector
        self.create_fake_file("parameters.md")
        self.create_fake_file("species.md")
        self.create_fake_file("reactions.md")

    def tearDown(self):
        # Clean up test directory
        if os.path.exists(self.test_directory):
            for root, dirs, files in os.walk(self.test_directory, topdown=False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
            os.rmdir(self.test_directory)

    def create_fake_file(self, file_name, content=None):
        """
        Helper function to create fake files with optional content.
        """
        if content is None:
            content = "# Fake content for {}\n".format(file_name)
        with open(os.path.join(self.test_directory, file_name), "w") as file:
            file.write(content)

    @patch("melektrodica.collector.DataParameters", autospec=True)
    @patch("melektrodica.collector.DataSpecies", autospec=True)
    @patch("melektrodica.collector.DataReactions", autospec=True)
    def test_initialization(
        self, MockDataReactions, MockDataSpecies, MockDataParameters
    ):
        # Mock the dependencies of the Collector class
        mock_parameters = MockDataParameters.return_value
        mock_species = MockDataSpecies.return_value
        mock_reactions = MockDataReactions.return_value

        # Create an instance of the Collector class
        collector = Collector(self.test_directory, self.writer)

        # Assert that the directory is correctly set
        self.assertEqual(collector.directory, self.test_directory)

        # Assert that the mocked parameters, species, and reactions are set
        self.assertIs(collector.parameters, mock_parameters)
        self.assertIs(collector.species, mock_species)
        self.assertIs(collector.reactions, mock_reactions)

        # Ensure mocks are called with the correct arguments
        MockDataParameters.assert_called_once_with(
            os.path.join(self.test_directory, "parameters.md"), self.writer
        )
        MockDataSpecies.assert_called_once_with(
            os.path.join(self.test_directory, "species.md"),
            mock_parameters,
            self.writer,
        )
        MockDataReactions.assert_called_once_with(
            os.path.join(self.test_directory, "reactions.md"),
            mock_parameters,
            mock_species,
            self.writer,
        )

    def test_raw_data_valid_file(self):
        # Create a valid test file
        test_file_path = os.path.join(self.test_directory, "test_file.md")
        with open(test_file_path, "w") as file:
            file.writelines(self.sample_file_content)

        header, data = Collector.raw_data(test_file_path, self.writer)
        self.assertEqual(header, ["Column1", "Column2", "Column3"])
        np.testing.assert_array_equal(data, [["1", "2", "3"], ["4", "5", "6"]])

    def test_raw_data_file_not_found(self):
        # Test when the file does not exist
        with self.assertRaises(FileNotFoundError):
            Collector.raw_data("non_existent_file.md", self.writer)

    def test_raw_data_empty_file(self):
        # Create an empty test file
        empty_file_path = os.path.join(self.test_directory, "empty_file.md")
        with open(empty_file_path, "w") as file:
            pass

        # Attempt to parse empty file
        with self.assertRaises(IndexError):
            Collector.raw_data(empty_file_path, self.writer)

    def test_column_exists_valid_column(self):
        # Test column exists functionality with a valid column
        Collector.column_exists("Column1", self.header, "test_file.md", self.writer)

    def test_column_exists_missing_column(self):
        # Test column exists functionality when column is missing
        with self.assertRaises(ValueError):
            Collector.column_exists(
                "MissingColumn", self.header, "test_file.md", self.writer
            )

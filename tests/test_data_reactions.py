import unittest
from unittest.mock import MagicMock, patch
import numpy as np
from melektrodica.collector import DataReactions


class TestDataReactionsMarkdown(unittest.TestCase):
    """
    Unit test class for DataReactions adapted to handle Markdown (.md) files
    containing reaction data.
    """

    class ParametersMock:
        def __init__(
            self,
            experimental=False,
            thermochemical=False,
            cstr=False,
            dg_reaction=False,
        ):
            self.experimental = experimental
            self.thermochemical = thermochemical
            self.cstr = cstr
            self.dg_reaction = dg_reaction

    class SpeciesMock:
        def __init__(self):
            self.list = ["H2", "O2", "H2O", "H+", "e-", "Pt", "O*", "OH*"]
            self.catalyst = []
            self.adsorbed = []

    def setUp(self):
        self.writer = MagicMock()
        self.parameters = self.ParametersMock()
        self.species = self.SpeciesMock()

    @staticmethod
    def parse_markdown_table(md_content):
        """
        Parses a Markdown table and returns its headers and data.

        :param md_content: A string containing the Markdown table.
        :return: A tuple (headers, rows) where headers is a list of column names
                 and rows is a numpy array with the table content.
        """
        # Use markdown-it or a custom approach to parse the table
        # Extract lines and differentiate headers and data
        lines = md_content.strip().split("\n")
        headers = [h.strip() for h in lines[0].split("|") if h.strip()]
        data = []
        for row in lines[2:]:  # Skip header and separator lines
            values = [v.strip() for v in row.split("|") if v.strip()]
            data.append(values)
        return headers, np.array(data)

    @patch("melektrodica.collector.Collector.column_exists")
    @patch("melektrodica.collector.Collector.raw_data")
    def test_initialization_with_markdown_file(self, mock_raw_data, mock_column_exists):
        # Simulate Markdown table content similar to the one you provided
        md_content = """
        | id | Reactions                    |    Ga | Beta |
        |----|------------------------------|------:|-----:|
        | DA | 0.5O2 + Pt <->  O*           | 0.391 |  0.0 |
        | RA | 0.5O2 + Pt + H+ + e- <-> OH* | 0.609 |  0.5 |
        | RT | O* + H+ + e- <-> OH*         | 0.590 |  0.5 |
        | RD | OH* + H+ + e- <-> H2O + Pt   | 0.278 |  0.5 |
        """

        # Mock the reading of raw Markdown data
        headers, rows = self.parse_markdown_table(md_content)
        mock_raw_data.return_value = (headers, rows)

        # Test initialization with parsed Markdown data
        self.parameters.thermochemical = True
        DataReactions(
            "test_reaction_file.md", self.parameters, self.species, self.writer
        )

        mock_column_exists.assert_any_call(
            "Ga", headers, "test_reaction_file.md", self.writer
        )
        self.writer.message.assert_called_with(
            "Thermochemical reactions parameters processed."
        )

    def test_process_reaction_valid(self):
        side = "0.5O2 + Pt"
        species_list = self.species.list
        result = DataReactions.process_reaction(side, species_list)
        self.assertEqual(result[0], ["O2", "Pt"])
        self.assertEqual(result[1], [0.5, 1.0])

    def test_process_reaction_invalid_species(self):
        side = "0.5O2 + Unknown"
        species_list = self.species.list
        with self.assertRaises(ValueError):
            DataReactions.process_reaction(side, species_list)

    def test_process_reaction_invalid_string_format(self):
        side = "O2 +++ Pt"
        species_list = self.species.list
        with self.assertRaises(ValueError):
            DataReactions.process_reaction(side, species_list)


if __name__ == "__main__":
    unittest.main()

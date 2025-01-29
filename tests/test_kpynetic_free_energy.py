"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Kpynetic, Unit test

"""

import unittest
import numpy as np
from unittest.mock import MagicMock
from melektrodica import *


class TestFreeEnergy(unittest.TestCase):
    """
    Unit test class for testing the FreeEnergy class.

    This class contains unit tests to verify the correctness of the FreeEnergy class
    by testing the 'reaction' method with various inputs. The tests include scenarios
    with positive, negative, and zero values in the upsilon and G_formation arrays.
    Additionally, it checks for the expected exceptions when the input arrays have
    mismatched lengths.

    Methods
    -------
    setUp():
        Initializes the FreeEnergy instance before each test.

    test_reaction_positive_values():
        Tests the reaction method with positive and negative values.

    test_reaction_negative_values():
        Tests the reaction method with all negative and positive values.

    test_reaction_zeros():
        Tests the reaction method with zero values.

    test_reaction_mismatched_lengths():
        Tests the reaction method with arrays of mismatched lengths.
    """

    def setUp(self):
        self.data = None
        self.FreeEnergy = kpynetic.FreeEnergy(self.data)

    def test_reaction_positive_values(self):
        upsilon = np.array([1, -1, -1])
        G_formation = np.array([-237.13, 0, -285.83])
        expected_result = -(upsilon @ G_formation)
        result = self.FreeEnergy.reaction(upsilon, G_formation)
        self.assertAlmostEqual(result, expected_result, places=2)

    def test_reaction_negative_values(self):
        upsilon = np.array([-1, 2, 1])
        G_formation = np.array([-237.13, 0, -285.83])
        expected_result = -(upsilon @ G_formation)
        result = self.FreeEnergy.reaction(upsilon, G_formation)
        self.assertAlmostEqual(result, expected_result, places=2)

    def test_reaction_zeros(self):
        upsilon = np.array([0, 0, 0])
        G_formation = np.array([-237.13, 0, -285.83])
        expected_result = 0.0
        result = self.FreeEnergy.reaction(upsilon, G_formation)
        self.assertAlmostEqual(result, expected_result, places=2)

    def test_reaction_mismatched_lengths(self):
        upsilon = np.array([1, -1])
        G_formation = np.array([-237.13, 0, -285.83])
        with self.assertRaises(ValueError):
            self.FreeEnergy.reaction(upsilon, G_formation)


if __name__ == "__main__":
    unittest.main()

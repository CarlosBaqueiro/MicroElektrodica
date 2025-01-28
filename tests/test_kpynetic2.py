import numpy as np
import unittest
from melektrodica.kpynetic import ReactionRate
from unittest.mock import MagicMock


class TestReactionRate(unittest.TestCase):
    class MockSpecies:
        def __init__(self):
            self.ns_catalyst = np.array([[0.5, 0.3]])

    class MockData:
        def __init__(self):
            self.species = TestReactionRate.MockSpecies()

    def setUp(self):
        self.mock_data = self.MockData()
        self.reaction_rate = ReactionRate(self.mock_data)

    def test_rate(self):
        k_rate = np.array([1.0, 2.0])
        c_reactants = np.array([[1.0, 2.0], [1.5, 2.5]])
        c_products = np.array([[0.5, 1.0], [0.8, 1.2]])
        theta = np.array([0.4, 0.6])
        upsilon = np.array([[2.0, -1.0], [-0.5, 3.0]])
        self.reaction_rate.concentrate = MagicMock(return_value=np.array([[1, 2], [3, 4]]))
        self.reaction_rate.power_law = MagicMock(return_value=np.array([[1, 5], [2, 6]]))

        result = self.reaction_rate.rate(k_rate, c_reactants, c_products, theta, upsilon)

        self.reaction_rate.concentrate.assert_called_once_with(c_reactants, c_products, theta)
        self.reaction_rate.power_law.assert_called_once_with(self.reaction_rate.concentrate.return_value, upsilon)
        np.testing.assert_array_equal(result, np.array([1 * 1 + 2 * 2, 5 * 1 + 6 * 2]))

    def test_empty_sites(self):
        theta = np.array([0.1, 0.2])
        expected_result = np.array([1 - 0.5 * 0.1 - 0.3 * 0.2])
        result = self.reaction_rate.empty_sites(theta)
        np.testing.assert_array_equal(result, expected_result)

    def test_concentrate(self):
        c_reactants = np.array([1.0, 2.0])
        c_products = np.array([0.5, 1.0])
        theta = np.array([0.4, 0.6])
        self.reaction_rate.empty_sites = MagicMock(return_value=np.array([0.7]))
        result = self.reaction_rate.concentrate(c_reactants, c_products, theta)
        expected_result = np.array([1.0, 2.0, 0.5, 1.0, 0.4, 0.6, 0.7])
        np.testing.assert_array_equal(result, expected_result)

    def test_power_law(self):
        concentration = np.array([[1.0, 2.0], [2.0, 3.0]])
        upsilon = np.array([[-1.0, 2.0], [1.0, -2.0]])
        expected_result = np.array([
            [np.prod([1.0 ** 1.0, 2.0 ** -2.0]), np.prod([2.0 ** -1.0, 3.0 ** 2.0])],
            [-np.prod([1.0 ** -1.0, 2.0 ** 2.0]), -np.prod([2.0 ** 1.0, 3.0 ** -2.0])]
        ])
        result = self.reaction_rate.power_law(concentration, upsilon)
        np.testing.assert_array_almost_equal(result, expected_result, decimal=6)


if __name__ == "__main__":
    unittest.main()

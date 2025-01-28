import numpy as np
import unittest
from melektrodica.kpynetic import ReactionRate
from unittest.mock import MagicMock


class TestReactionRate(unittest.TestCase):
    class MockSpecies:
        def __init__(self):
            self.ns_catalyst = np.array([[1, 2]])

    class MockData:
        def __init__(self):
            self.species = TestReactionRate.MockSpecies()

    def setUp(self):
        self.mock_data = self.MockData()
        self.reaction_rate = ReactionRate(self.mock_data)

    def test_rate(self):
        k_rate = np.array([1.0, 1000.0])
        c_reactants = np.array([1.5, 2.5])
        c_products = np.array([0.8, 1.2])
        theta = np.array([0.2, 0.3])
        upsilon = np.array([-1.0, -2.0, 1.0, 0.0, 2.0, 2.0, -1.0])
        empty_sites = np.array([1-1.0*0.2-2.0*0.3])
        concentration = np.concatenate((c_reactants, c_products, theta, empty_sites))
        self.reaction_rate.concentrate = MagicMock(return_value=concentration)
        power_law = np.array([1.5**1*2.5**2*empty_sites[0]**1, 0.8**1*1.2**0*0.2**2*0.3**2])
        self.reaction_rate.power_law = MagicMock(return_value=power_law)
        result = self.reaction_rate.rate(k_rate, c_reactants, c_products, theta, upsilon)
        self.reaction_rate.concentrate.assert_called_once_with(c_reactants, c_products, theta)
        self.reaction_rate.power_law.assert_called_once_with(self.reaction_rate.concentrate.return_value, upsilon)
        np.testing.assert_array_equal(result, np.sum(k_rate * power_law))

    def test_powerlaw_positive_upsilon(self):
        concentrations = np.array([[1.0, 2.0], [2.0, 3.0]])
        upsilon = np.array([1.0, 1.0])
        result = self.reaction_rate.power_law(concentrations, upsilon)
        self.assertTrue(np.all(result[1] < 0))

    def test_powerlaw_negative_upsilon(self):
        concentrations = np.array([[1.0, 2.0], [2.0, 3.0]])
        upsilon = np.array([-1.0, -1.0])
        result = self.reaction_rate.power_law(concentrations, upsilon)
        self.assertTrue(np.all(result[0] > 0))

    def test_empty_sites(self):
        theta = np.array([0.1, 0.2])
        expected_result = np.array([1 - 1.0 * 0.1 - 2.0 * 0.2])
        result = self.reaction_rate.empty_sites(theta)
        np.testing.assert_array_equal(result, expected_result)

    def test_concentrate_combines_inputs(self):
        c_reactants = np.array([1.0, 2.0])
        c_products = np.array([0.5, 0.5])
        theta = np.array([0.1, 0.2])
        result = self.reaction_rate.concentrate(c_reactants, c_products, theta)
        self.assertEqual(
            result.shape[0], c_reactants.size + c_products.size + theta.size + 1
        )

    def test_concentrate(self):
        c_reactants = np.array([1.0, 2.0])
        c_products = np.array([0.5, 1.0])
        theta = np.array([0.4, 0.6])
        self.reaction_rate.empty_sites = MagicMock(return_value=np.array([0.7]))
        result = self.reaction_rate.concentrate(c_reactants, c_products, theta)
        expected_result = np.array([1.0, 2.0, 0.5, 1.0, 0.4, 0.6, 0.7])
        np.testing.assert_array_equal(result, expected_result)

    def test_power_law(self):
        concentration = np.array([0.1, 2.0])
        upsilon = np.array([[-1.0, 2.0],
                            [3.0, -4.0]])
        expected_result = np.array([
            [np.prod([0.1 ** 1.0]), -np.prod([2.0 ** 2.0])],
            [np.prod([2.0 ** 4.0]), -np.prod([0.1 ** 3.0])]
        ]).T
        result = self.reaction_rate.power_law(concentration, upsilon)
        np.testing.assert_array_almost_equal(result, expected_result, decimal=6)


if __name__ == "__main__":
    unittest.main()

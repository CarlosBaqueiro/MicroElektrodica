"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Kpynetic, Unit test

"""

import unittest
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
        self.FreeEnergy = FreeEnergy()

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


class TestRateConstants(unittest.TestCase):
    """
    Test cases for the RateConstants class, verifying the correctness of its methods.

    This unittest class includes tests for different methods in the RateConstants
    class: constant, experimental, chemical, and electrical. Each test case checks
    the accuracy of the computations and ensures that the methods return the expected
    results under specific conditions.
    """

    def setUp(self):
        self.rate_constants = RateConstants()

    def test_constant(self):
        class Parameters:
            T = 300

        self.rate_constants.data = self
        self.rate_constants.data.parameters = Parameters()

        pre_exponential = 2
        experimental = 2
        chemical = 10
        electrical = 5
        expected_result = (
            pre_exponential
            * experimental
            * np.exp(-(chemical + electrical) / (8.617333262145e-5) / 300)
        )

        result = self.rate_constants.constant(
            pre_exponential, experimental, chemical, electrical
        )
        self.assertAlmostEqual(result, expected_result, places=6)

    def test_experimental(self):
        forward_constants = 5
        backward_constants = 3
        expected_result = np.array([forward_constants, backward_constants])
        result = RateConstants.experimental(forward_constants, backward_constants)
        np.testing.assert_array_equal(result, expected_result)

    def test_chemical(self):
        G_activation = 100
        DG_reaction = 80
        ne = 1
        expected_result = np.array([G_activation, G_activation + DG_reaction])
        result = RateConstants.chemical(G_activation, DG_reaction, ne)
        np.testing.assert_array_equal(result, expected_result)

    def test_electrical(self):
        eta = 0.1
        ne = 2
        beta = 0.5
        expected_result = np.array([-ne * beta * eta, ne * (1 - beta) * eta])
        result = RateConstants.electrical(eta, ne, beta)
        np.testing.assert_array_equal(result, expected_result)


class TestReactionRate(unittest.TestCase):
    """
    Test suite for the ReactionRate class.

    This class contains unit tests for the ReactionRate class to ensure its methods
    function correctly. It uses the unittest framework and creates a ReactionRate
    instance before each test to check various aspects of the class methods.

    Methods:
        setUp: Initializes the test setup including creating a ReactionRate instance.
        test_emptysites_valid_theta: Checks if the emptysites method returns values
            all less than 1.0 for a given theta.
        test_concentrate_combines_inputs: Verifies that the concentrate method
            combines input arrays correctly and checks the shape of the result.
        test_powerlaw_positive_upsilon: Tests the powerlaw method when upsilon is
            positive to ensure the results meet expected conditions.
        test_powerlaw_negative_upsilon: Tests the powerlaw method when upsilon is
            negative to ensure the results meet expected conditions.
    """

    def setUp(self):
        self.reaction_rate = ReactionRate()

        class Species:
            nu_catalyst = np.array([[0.5], [0.2]])

        self.reaction_rate.species = Species()

    def test_emptysites_valid_theta(self):
        theta = np.array([0.2])
        result = self.reaction_rate.emptysites(theta)
        self.assertTrue((result < 1.0).all())

    def test_concentrate_combines_inputs(self):
        c_reactants = np.array([1.0, 2.0])
        c_products = np.array([0.5, 0.5])
        theta = np.array([0.1])
        result = self.reaction_rate.concentrate(c_reactants, c_products, theta)
        self.assertEqual(
            result.shape[0], c_reactants.size + c_products.size + theta.size + 2
        )

    def test_powerlaw_positive_upsilon(self):
        concentrations = np.array([[1.0, 2.0], [2.0, 3.0]])
        upsilon = np.array([1.0, 1.0])
        result = self.reaction_rate.powerlaw(concentrations, upsilon)
        self.assertTrue(np.all(result[1] < 0))

    def test_powerlaw_negative_upsilon(self):
        concentrations = np.array([[1.0, 2.0], [2.0, 3.0]])
        upsilon = np.array([-1.0, -1.0])
        result = self.reaction_rate.powerlaw(concentrations, upsilon)
        self.assertTrue(np.all(result[0] > 0))


class TestKpynetic(unittest.TestCase):
    """
    A unit test class for the Kpynetic module.

    This class contains unit tests that verify the correct initialization and behavior of the
    Kpynetic module, ensuring that attributes and methods function as expected.
    """

    class Parameters:
        anode = False
        pre_exponential = 1.0
        js = False
        js_value = 1.0
        tst = False
        kappa = 1.0
        T = 300
        m = 1.0
        experimental = False
        chemical = False
        DG_reaction = False
        G_formation = False
        potential = 0.0

    class Species:
        nu_catalyst = np.array([[0.5], [0.2]])

    class Reactions:
        list = []
        nua = np.array([1.0])
        ne = np.array([1.0])
        beta = np.array([1.0])
        nu = np.array([-1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0])

    class Data:
        def __init__(self):
            self.parameters = TestKpynetic.Parameters()
            self.species = TestKpynetic.Species()
            self.reactions = TestKpynetic.Reactions()

    def setUp(self):
        self.data = TestKpynetic.Data()
        self.kpynetic = Kpynetic(self.data)

    def test_initialization(self):
        self.assertEqual(self.kpynetic.electrode, -1.0)
        self.assertEqual(self.kpynetic.pre_exp, 1.0)
        np.testing.assert_array_equal(
            self.kpynetic.eK, np.ones((2, len(self.data.reactions.list)))
        )
        np.testing.assert_array_equal(
            self.kpynetic.qK, np.zeros((2, len(self.data.reactions.list)))
        )


if __name__ == "__main__":
    unittest.main()

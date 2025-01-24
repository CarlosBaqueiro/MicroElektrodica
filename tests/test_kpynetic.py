# tests/test_kpynetic.py

import numpy as np
import unittest
from melektrodica.kpynetic import Kpynetic


class TestKpynetic(unittest.TestCase):
    """
    Unit test class for Kpynetic to validate the functionality of its methods
    including normal use cases, edge cases, and exception scenarios.
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
        thermochemical = False
        DG_reaction = False
        G_formation = False
        potential = 0.0

    class Species:
        nu_catalyst = np.array([[0.5], [0.2]])

    class Reactions:
        list = []
        upsilon_a = np.array([1.0])
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

    def test_foverpotential_valid(self):
        """Test that the foverpotential method produces correct outputs for valid inputs."""
        potential = 0.75
        c_reactants = np.array([1.0, 1.5])
        c_products = np.array([0.8, 1.2])
        theta = np.array([0.1, 0.05])
        result = self.kpynetic.foverpotential(potential, c_reactants, c_products, theta)
        self.assertIsInstance(result, float)

    def test_foverpotential_mismatched_array_sizes(self):
        """Ensure that foverpotential raises a ValueError for mismatched array sizes."""
        potential = 0.75
        c_reactants = np.array([1.0])
        c_products = np.array([1.2, 0.8])
        theta = np.array([0.1, 0.05])
        with self.assertRaises(ValueError):
            self.kpynetic.foverpotential(potential, c_reactants, c_products, theta)

    def test_get_argument_valid(self):
        """Test that get_argument computes correctly with a valid potential."""
        potential = 1.2
        argument = self.kpynetic.get_argument(potential)
        self.assertIsInstance(argument, float)

    def test_get_argument_invalid(self):
        """Ensure get_argument raises an exception for invalid input."""
        potential = None  # An invalid input
        with self.assertRaises(TypeError):
            self.kpynetic.get_argument(potential)

    def test_current_valid(self):
        """Test that the current method works correctly with valid inputs."""
        potential = 0.5
        c_reactants = np.array([1.0, 0.8])
        c_products = np.array([0.6, 0.4])
        theta = np.array([0.1, 0.2])
        result = self.kpynetic.current(potential, c_reactants, c_products, theta)
        self.assertIsInstance(result, np.ndarray)

    def test_current_invalid_array_shape(self):
        """Ensure that current raises a ValueError for invalid array shapes."""
        potential = 0.5
        c_reactants = np.array([1.0])
        c_products = np.array([0.6, 0.4])
        theta = np.array([0.1, 0.2])
        with self.assertRaises(ValueError):
            self.kpynetic.current(potential, c_reactants, c_products, theta)

    def test_dcdt_valid(self):
        """Test correct computation of dcdt with valid inputs."""
        rate = np.array([1.0, -1.0, 0.5])
        upsilon = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
        result = self.kpynetic.dcdt(rate, upsilon)
        self.assertIsInstance(result, np.ndarray)

    def test_dcdt_invalid_upsilon_shape(self):
        """Ensure that dcdt raises an exception for invalid upsilon array shape."""
        rate = np.array([1.0, -1.0])
        upsilon = np.array([[1, -1]])  # Invalid shape for multiplication
        with self.assertRaises(ValueError):
            self.kpynetic.dcdt(rate, upsilon)

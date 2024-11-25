"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Kpynetic class

"""

import numpy as np
from numpy import ndarray

from .Constants import F, k_B, h

# for debugging
#from .Tools import showme
#import sys
# sys.exit()

class FreeEnergy:
    """
    Provides static methods for calculating different types of free energy changes.

    Detailed purpose includes calculations like reaction free energy, reduction reaction free energy,
    and activation reaction free energy typically used in thermodynamics and electrochemistry contexts.

    """

    def __init__(self):
        pass

    @staticmethod
    def reaction(upsilon: ndarray, G_formation: ndarray) -> ndarray:
        """
        Calculate the Gibbs free energy change for a reaction.

        This method computes the negative of the dot product between the
        vector of stoichiometric coefficients and the vector of Gibbs free
        energies of formation for each substance involved in the reaction.

        :param upsilon: A vector of stoichiometric coefficients
                        representing the number of moles of each substance
                        participating in the reaction.
        :type upsilon: numpy.ndarray

        :param G_formation: A vector of Gibbs free energies of formation
                            for each substance involved in the reaction.
        :type G_formation: numpy.ndarray

        :return: The Gibbs free energy change for the reaction.
        :rtype: float
        """
        return -(upsilon @ G_formation)


class RateConstants:
    """
    A class used to calculate various rate constants for a chemical reaction system.

    This class provides methods to compute general rate constants, experimental rate constants,
    chemical reaction factors, and electrical parameters.

    :ivar parameters: Parameters object containing temperature and other relevant properties.
    :type parameters: object
    """

    def __init__(self):
        pass

    def constant(self, pre_exponential=1, experimental=1, chemical=0, electrical=0):
        """
        Example class for demonstrating the calculation of a kinetics constant based on
        provided parameters.
        :param pre_exponential: Preexponential factor.
        :param experimental: Experimental rate constant.
        :param chemical: Chemical part of the rate constant.
        :param electrical: Electrical part of the rate constant.
        :param parameters: Object with temperature attribute
        :type parameters: Object
        """
        argument = (chemical + electrical)
        return pre_exponential * experimental * np.exp(-argument / k_B / self.data.parameters.T)

    @staticmethod
    def experimental(forward_constants, backward_constants):
        """
        Converts the forward and backward constants into a numpy array.

        This function takes the provided forward and backward constants,
        and combines them into a single numpy array. This is useful
        for numerical computations where a unified data structure
        for constants is needed.

        :param forward_constants: Constants used in the forward computation
        :type forward_constants: Any
        :param backward_constants: Constants used in the backward computation
        :type backward_constants: Any
        :return: A numpy array combining forward and backward constants
        :rtype: numpy.ndarray
        """
        return np.array([forward_constants, backward_constants])

    @staticmethod
    def chemical(G_activation, DG_reaction, ne):
        """
        Compute the chemical energy.

        This method calculates the chemical energy based on the activation
        energy, change in the Gibbs free energy of the reaction, and the
        number of electrons involved.

        :param G_activation: Activation energy.
        :param DG_reaction: Change in Gibbs free energy of the reaction.
        :param ne: Number of electrons involved.
        :return: Array containing the activation energy and the total chemical
                 energy (activation energy plus change in Gibbs free energy).
        """
        return np.array([G_activation, G_activation + DG_reaction])

    @staticmethod
    def electrical(eta, ne, beta):
        """
        Calculates the electrical field components based on the specified parameters.

        The function utilizes the given eta, ne, and beta values to compute an array
        that represents the electrical field. The resulting array contains two elements:
        the first being the product of -ne, beta, and eta, and the second being the product
        of ne, (1 - beta), and eta.

        :param eta: The efficiency factor.
        :type eta: float
        :param ne: The charge density.
        :type ne: float
        :param beta: The proportionality constant.
        :type beta: float
        :return: An array representing the electrical field components.
        :rtype: numpy.ndarray
        """
        return np.array([-ne * beta * eta, ne * (1 - beta) * eta])


class ReactionRate:
    """
    Computes the rate of a reaction based on various parameters.

    The `ReactionRate` class calculates the reaction rate using the provided rate constants,
    reactant and product concentrations, surface coverage, and reaction orders.

    :ivar species: Species involved in the reaction.
    :type species: SomeSpeciesClass

    :ivar nu_catalyst: Stoichiometric coefficients for the catalyst.
    :type nu_catalyst: np.ndarray
    """

    def __init__(self):
        pass

    def rate(self, k_rate, c_reactants, c_products, theta, upsilon):
        """
        Calculate the reaction rate using given parameters.

        This method computes the rate of reaction based on the provided
        rate constant, concentrations of reactants and products, a
        parameter theta, and a power law exponent upsilon.

        :param k_rate: The rate constant for the reaction.
        :param c_reactants: The concentrations of the reactants.
        :param c_products: The concentrations of the products.
        :param theta: A parameter used to calculate concentrations.
        :param upsilon: The exponent used in the power law for rate
                        computation.
        :return: The calculated reaction rate.
        """
        concentrations = self.concentrate(c_reactants, c_products, theta)
        rate = self.powerlaw(concentrations, upsilon)
        return np.sum(k_rate * rate, axis=0)

    def emptysites(self, theta):
        """
        Calculate and return the empty sites in the catalyst.

        This function computes the empty sites in the catalyst by performing a matrix
        multiplication between the catalyst species matrix (`nu_catalyst`) and the
        input array `theta`, and then subtracting the result from 1.

        :param theta: Array representing the site occupancy mentioned in the
          catalyst species matrix.
        :type theta: np.ndarray
        :return: Array representing the empty sites in the catalyst.
        :rtype: np.ndarray
        """
        return np.array(1 - self.species.nu_catalyst @ theta)

    def concentrate(self, c_reactants, c_products, theta):
        """
        This method concatenates the concentrations of reactants, products, the given theta value,
        and the result of `self.emptysites()` method into a single numpy array.

        :param c_reactants: Concentration vector for reactants.
        :type c_reactants: np.ndarray
        :param c_products: Concentration vector for products.
        :type c_products: np.ndarray
        :param theta: Theta value used for calculations.
        :type theta: np.ndarray
        :return: A concatenated numpy array containing all input concentrations and the result of
                 the `self.emptysites()` method.
        :rtype: np.ndarray
        """
        return np.concatenate([c_reactants, c_products, theta, self.emptysites(theta)])

    def powerlaw(self, concentration: ndarray, upsilon: ndarray) -> ndarray:
        """
        Computes the power law of given concentrations with a specified upsilon value.

        This method calculates two arrays based on whether upsilon is negative or positive, raising
        the concentration array to the power of upsilon and taking the product across the specified axis.

        :return:
        :param concentration: Array of concentrations on which power law is computed
        :type concentration: np.ndarray
        :param upsilon: Exponent value used in the power law calculation
        :type upsilon: float

        :return: Array containing the computed power law values
        :rtype: np.ndarray
        """
        return np.array([np.prod(concentration ** (-upsilon * (upsilon < 0)), axis=1),
                         -np.prod(concentration ** (upsilon * (upsilon > 0)), axis=1)])


class Kpynetic(FreeEnergy, RateConstants, ReactionRate):
    """
    Kpynetic class for modeling and simulating kinetics of chemical reactions.

    This class integrates methods from FreeEnergy, RateConstants, and ReactionRate
    to calculate reaction rates, current, and concentration changes over time. It
    is designed to handle various parameters, including pre-exponential factors,
    experimental rate constants, and chemical rate constants.

    :ivar data: Input data containing parameters, species, and reactions.
    :type data: Data
    :ivar parameters: Parameters for the reactions.
    :type parameters: Parameters
    :ivar species: Species involved in the reactions.
    :type species: Species
    :ivar reactions: Reaction information.
    :type reactions: Reactions
    :ivar electrode: Electrode potential sign, set based on parameters.
    :type electrode: float
    :ivar pre_exp: Pre-exponential factor, calculated based on transition state
        theory (TST) or predefined value.
    :type pre_exp: float
    :ivar eK: Experimental rate constants array.
    :type eK: np.ndarray
    :ivar qK: Chemical rate constants array.
    :type qK: np.ndarray
    :ivar Ga: Gibbs free energy of activation, set if parameters.chemical is True.
    :type Ga: np.ndarray or None
    :ivar DG_reaction: Gibbs free energy change of reaction, set if parameters.chemical
        is True.
    :type DG_reaction: np.ndarray or None
    """

    def __init__(self, data):
        self.data = data
        self.parameters = data.parameters
        self.species = data.species
        self.reactions = data.reactions

        self.electrode = 1.0
        if not self.parameters.anode: self.electrode = - 1.0

        # Pre-exponential
        self.pre_exp = self.parameters.pre_exponential
        if self.parameters.js:
            self.pre_exp = self.parameters.js_value / F
        if self.parameters.tst:
            self.pre_exp = self.parameters.kappa * k_B * self.parameters.T ** self.parameters.m / h

        # Experimental
        self.eK = np.ones((2, len(self.reactions.list)))
        if self.parameters.experimental:
            self.eK = RateConstants.experimental(self.reactions.k_f, self.reactions.k_b)

        # Chemical
        self.qK = np.zeros((2, len(self.reactions.list)))
        if self.parameters.chemical:
            self.Ga = self.reactions.Ga
            if self.parameters.DG_reaction:
                self.DG_reaction = self.reactions.DG_reaction
            elif self.parameters.G_formation:
                self.G_formation = self.species.G_formation_ads
                self.DG_reaction = FreeEnergy.reaction(self.reactions.nua, self.G_formation)
            self.qK = RateConstants.chemical(self.Ga, self.DG_reaction, self.reactions.ne)

    def potential_function(self, potential, c_reactants, c_products, theta):
        eta = potential
        self.potentialK = self.electrode * RateConstants.electrical(eta, self.reactions.ne, self.reactions.beta)
        self.k_rate = self.constant(pre_exponential=self.pre_exp,
                                    experimental=self.eK, chemical=self.qK,
                                    electrical=self.potentialK)
        self.v = self.rate(self.k_rate, c_reactants, c_products, theta, self.reactions.nu)

    def current(self, potential, c_reactants, c_products, theta):
        self.potential_function(potential, c_reactants, c_products, theta)
        return np.dot(self.reactions.ne, self.v) * F

    def dcdt(self, rate, upsilon: np.ndarray) -> np.ndarray:
        return np.dot(rate, upsilon)

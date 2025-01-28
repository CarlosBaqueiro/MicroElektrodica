"""

    μElektrodica© 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Kpynetic class

"""

import copy
import numpy as np
from numpy import ndarray

from .writer import Writer
from .constants import F, k_B, h


# for debugging
# from .Tools import showme
# import sys
# sys.exit()


class FreeEnergy:
    """
    Class representing an entity with initialization using the given data.

    The class accepts a single argument during initialization and sets up
    an internal representation of the provided data. This is typically
    used as a base or utility class for more complex operations or data
    management systems.

    Attributes:
        data: The data object provided during initialization, which is
            stored internally for further operations or usage.
    """

    def __init__(self, data: object):
        """
        Class representing an entity with initialization using the given data.

        The class accepts a single argument during initialization and sets up
        an internal representation of the provided data. This is typically
        used as a base or utility class for more complex operations or data
        management systems.

        Attributes:
            data: The data object provided during initialization, which is
                stored internally for further operations or usage.

        """
        pass

    @staticmethod
    def reaction(upsilon: ndarray, g_formation: ndarray) -> ndarray:
        """
        Calculate the Gibbs free energy change for a reaction.

        This method returns the reaction standard Gibbs free energy of the step reaction :math:`i` :cite:p:`SmithThermo`

        .. math::

            \\Delta G^\\circ_r = -\\sum_j \\upsilon_{ij} G_j^\\circ

        where :math:`G_j^\\circ` is the standard Gibbs free formation energy of :math:`j` species

        :param upsilon: A vector of stoichiometric coefficients
                        representing the number of moles of each substance
                        participating in the reaction.
        :type upsilon: numpy.ndarray

        :param g_formation: A vector of Gibbs free energies of formation
                            for each substance involved in the reaction.
        :type g_formation: numpy.ndarray

        :return: The Gibbs free energy change for the reaction.
        :rtype: float
        """
        # print(f'Kpynetic.FreeEnergy.reaction: \n upsilon = {upsilon}, g_formation = {g_formation}')
        return -(upsilon @ g_formation)


class RateConstants:
    """
    Class for performing calculations related to chemical kinetics and thermodynamic properties.

    Kinetic rate constants, for an electrochemical system, can be determined experimentally or estimated
    by thermodynamic and electrochemical parts, the forward, :math:`\\overrightarrow{k_i}`,
    and backward, :math:`\\overleftarrow{k_i}`, rate constants :cite:p:`Exner2022`:

    .. math::

       \\overrightarrow{k_i} = A
            \\underbrace{\\exp\\left(\\frac{- \\Delta G^{\\circ}_{a,i}} {k_BT} \\right)}_{\\textit{Thermodynamic}}
            \\underbrace{\\exp\\left(\\frac{n_i\\beta_i\\eta} {k_BT} \\right)}_{\\textit{Electronic}}

    .. math::

        \\overleftarrow{k_i} =
            \\frac{\\overrightarrow{k_i}}{K_{eq,i}} = A
                \\underbrace{\\exp\\left(\\frac{-\\Delta G^{\\circ}_{a,i} +
                    \\Delta G^\\circ_{r,i}}{k_BT} \\right)}_{\\textit{Thermodynamic}}
                \\underbrace{\\exp\\left(\\frac{-n_i(1-\\beta_i)\\eta}{k_BT} \\right)}_{\\textit{Electronic}}


    This class is designed to facilitate the computation of various constants and parameters
    associated with chemical reactions. It uses data provided during initialization to access
    attributes like temperature, reactions, and other necessary parameters.

    Attributes
    ----------
    data : object
        Data object containing parameters and reactions.
    parameters : object
        Parameters object extracted from the input data.
    reactions : object
        Reaction information extracted from the input data.
    """

    def __init__(self, data):
        self.data = data
        self.parameters = data.parameters
        self.reactions = data.reactions

    def constant(
            self, pre_exponential=1, experimental=1, thermochemical=0, electronic=0
    ) -> ndarray:
        """
        Calculates a reaction rate constants using the provided parameters.

        This method computes the rate constant using the pre-exponential factor,
        experimental factor, thermochemical and electronic contributions.

        Parameters
        ----------
        pre_exponential : float, optional
            The pre-exponential factor in the Arrhenius equation, often denoted as A.
        experimental : float, optional
            A multiplicative factor representing experimental data adjustments.
        thermochemical : float, optional
            The thermochemical contribution to the activation energy.
        electronic : float, optional
            The electronic contribution to the activation energy.

        Returns
        -------
        float
            The calculated reaction rate constant.
        """
        self.argument = thermochemical + electronic
        return (
                pre_exponential
                * experimental
                * np.exp(-self.argument / k_B / self.data.parameters.temperature)
        )

    @staticmethod
    def experimental(
            forward_constants: np.ndarray, backward_constants: np.ndarray
    ) -> ndarray:
        """
        Construct an array combining forward and backward constants.

        The method compiles the given forward and backward constants into a
        NumPy array. This array retains the provided constants in their
        respective order for direct usage or subsequent computations.

        Parameters
        ----------
        forward_constants : ndarray
            Values to use as forward constants, of any data type appropriate
            for inclusion in an array.
        backward_constants : ndarray
            Values to use as backward constants, of any data type appropriate
            for inclusion in an array.

        Returns
        -------
        numpy.ndarray
            A NumPy array containing the provided forward constants and
            backward constants, in the given order.
        """
        return np.array([forward_constants, backward_constants])

    @staticmethod
    def thermochemical(
            g_activation: np.ndarray, g_formation: np.ndarray, upsilon: np.ndarray
    ) -> ndarray:
        """
        Calculates the thermochemical free energies based on activation energy, formation
        energy, and stoichiometric parameter.

        This static method computes the Gibbs free energy of the reaction using a given
        stoichiometric coefficient and formation energy. It then returns an array containing
        the activation energy and the sum of the activation energy and the Gibbs free energy
        of the reaction.

        Parameters
        ----------
        g_activation : float
            The activation energy of the reaction in appropriate energy units.
        g_formation : float
            The Gibbs free energy of formation in appropriate energy units.
        upsilon : float
            The stoichiometric coefficient for the reaction, a dimensionless parameter.

        Returns
        -------
        numpy.ndarray
            A 1D array where:
                - The first element is the activation energy (`g_activation`).
                - The second element is the sum of the activation energy and the Gibbs free
                  energy of the reaction.

        """

        dg_reaction = FreeEnergy.reaction(upsilon, g_formation)
        return np.array([g_activation, g_activation + dg_reaction])

    @staticmethod
    def electronic(eta: float, ne: np.ndarray, beta: np.ndarray) -> ndarray:
        """
        Computes the electronic energy and its complementary component based on the input parameters.

        This function calculates the electronic energy and its complementary component
        using the provided eta, ne, and beta parameters. It combines the number of
        electrons (ne), a scaling factor (beta), and an energy term (eta) to compute
        these values. The resulting array includes negative energy and complementary
        part in specific proportions.

        Parameters
        ----------
        eta : float
            Energy term used in the computation.
        ne : int
            Number of electrons involved in the computation.
        beta : float
            Scaling factor determining proportions of energy allocation.

        Returns
        -------
        numpy.ndarray
            A 1D array where the first element is the negative energy component,
            and the second element is the complementary energy component.
        """
        return np.array([-ne * beta * eta, ne * (1 - beta) * eta])


class ReactionRate:
    """
    Class to model reaction kinetics and catalysis.

    This class provides functionalities to calculate reaction rates based on
    the law of mass action and power law kinetics. It is designed to handle
    reaction rate constants, reactants and product concentrations, and
    catalyst site occupancy to compute reaction mechanisms effectively :cite:p:`motagamwala2018microkinetic,
    Motagamwala2020, bieniasz2015, razdan2023`.

    .. math::

          \\nu_i = \\overrightarrow{k_i}
                    \\prod_{\\substack{j \\\\ \\upsilon_{ij}<0}} c_j ^{-\\upsilon_{ij}}
                    \\prod_{\\substack{j \\\\ \\upsilon_{ij}<0}} \\theta_j ^{-\\upsilon_{ij}}
                    -
                    \\overleftarrow{k_i}
                    \\prod_{\substack{j \\\\ \\upsilon_{ij}>0}} c_j ^{\\upsilon_{ij}}
                    \\prod_{\substack{j \\\\ \\upsilon_{ij}>0}} \\theta_j ^{\\upsilon_{ij}}

    Attributes
    ----------
    data : object
        A data source containing information about the chemical system.
    species : object
        A property of the given data that includes catalyst species matrix or
        related information.
    """

    def __init__(self, data: object):
        self.data = data
        self.species = data.species

    def rate(self, k_rate, c_reactants, c_products, theta, upsilon) -> ndarray:
        """
        Calculate the reaction rate for a given reaction system.

        This method calculates the rate of reaction using the power-law
        kinetics model. The reaction rate is determined by computing the
        concentration of reactants and products based on system parameters,
        applying the power-law formula and multiplying by the rate constant.

        Parameters
        ----------
        k_rate : np.ndarray
            The reaction rate constants for the chemical reactions.
        c_reactants : np.ndarray
            The concentration of the reactants involved in the reactions.
        c_products : np.ndarray
            The concentration of the products resulting from the reactions.
        theta : np.ndarray
            The set of system parameters governing reactant and product behaviors.
        upsilon : np.ndarray
            Stoichiometric coefficients dictating reactant-product relationships.

        Returns
        -------
        np.ndarray
            The computed reaction rate for each reaction across the system.
        """

        concentrations = self.concentrate(c_reactants, c_products, theta)
        rate = self.power_law(concentrations, upsilon)
        return np.sum(k_rate * rate, axis=0)

    def empty_sites(self, theta: np.ndarray) -> ndarray:
        """
        Computes the empty sites on a catalytic surface.

        The function calculates the empty sites on a catalytic surface based on the
        provided theta vector, which represents the coverage of various species on the
        surface. The species' contribution is determined by the ns_catalyst matrix,
        which indicates the stoichiometric coefficients corresponding to the catalyst
        interaction.

        Parameters
        ----------
        theta : numpy.ndarray
            A 1-dimensional array representing the coverage of the species on the
            catalytic surface.

        Returns
        -------
        numpy.ndarray
            A 1-dimensional array showing the fraction of empty sites on the catalytic
            surface.
        """

        return np.array(1 - self.species.ns_catalyst @ theta)

    def concentrate(
            self, c_reactants: np.ndarray, c_products: np.ndarray, theta: np.ndarray
    ) -> ndarray:
        """
        Combine the concentrations of reactants, products, and other parameters into a single array.

        This method takes the concentrations of reactants and products, as well as a given
        parameter `theta`, and combines them with the result from the `empty_sites` method.
        The final output is a concatenated array that includes all these components.

        Parameters
        ----------
        c_reactants : np.ndarray
            Concentration array for reactants.
        c_products : np.ndarray
            Concentration array for products.
        theta : np.ndarray
            Array representing coverage or occupancy values.

        Returns
        -------
        np.ndarray
            A concatenated array containing the concentrations of reactants, products, the
            provided `theta`, and the result of the `empty_sites` method.
        """

        return np.concatenate([c_reactants, c_products, theta, self.empty_sites(theta)])

    def power_law(self, concentration: ndarray, upsilon: ndarray) -> ndarray:
        """
        Calculates the power law product for two given parameters (concentration and upsilon).

        The method computes an array with two rows. The first row is the product of the
        concentration raised to the power of the negative upsilon, conditioned by upsilon
        being negative. The second row is the product of the concentration raised to the
        power of the upsilon, conditioned by upsilon being positive.

        Parameters
        ----------
        concentration : ndarray
            A 2D array where each row corresponds to a set of concentration values.
        upsilon : ndarray
            A 2D array where each row corresponds to a set of upsilon values.

        Returns
        -------
        ndarray
            A 2D array where the first row corresponds to the computed power law products
            of concentration and negative upsilon under specific conditions, and the second
            row corresponds to the computed power law products of concentration and positive
            upsilon under specific conditions.
        """

        return np.array(
            [
                np.prod(concentration ** (-upsilon * (upsilon < 0)), axis=1),
                -np.prod(concentration ** (upsilon * (upsilon > 0)), axis=1),
            ]
        )


class Kpynetic(FreeEnergy, RateConstants, ReactionRate):
    """
    Represents a computational model for chemical kinetics and electrochemical dynamics.

    This class is designed to calculate the chemical kinetics and electrochemical behavior
    of reactions by combining thermochemical, experimental, and electronic contributions.
    It utilizes rate constants, free energy, and reaction parameters to compute reaction
    rates, overpotentials, and currents. The purpose is to analyze and simulate reaction
    mechanisms in various chemical and electrochemical systems.

    Attributes
    ----------
    writer : Writer
        Object used for logging and managing output messages.
    data : copy.deepcopy
        Deep copy of the input data for initializing the model.
    parameters : object
        Reaction parameters extracted from the input data.
    species : object
        List of species involved in the reactions.
    reactions : object
        Reactions data including forward and backward reaction rates.
    operation : object
        Operational parameters derived from input data.
    electrode : float
        Factor to determine the nature of the electrode (e.g., anode or cathode).
    pre_exp : float
        Pre-exponential factor used in rate calculations.
    experimental_part : np.ndarray
        Contributions of experimental factors to the reaction rates.
    thermochemical_part : np.ndarray
        Contributions of thermochemical factors to the reaction rates.
    g_activation : list or None
        Activation free energy values for the reactions.
    g_formation : list or None
        Formation free energy values for the adsorbed species.
    dg_reaction : list
        Reaction free energy changes computed from stoichiometry.
    electronic_part : float or None
        Electronic contribution to rate constants.
    constant : callable
        Method for evaluating the overall rate constant.
    rate : callable
        Method for computing reaction rates.
    v : np.ndarray or None
        Reaction rates computed during runtime.
    """

    def __init__(self, data: object, writer=None):
        """
        Initializes the kinetic model object for reaction simulation based on given data.

        This class is responsible for initializing and setting up various reaction-based
        parameters and constants required for kinetic simulations. The parameters include
        pre-exponential factors, electrode properties, experimental reaction constants, and
        thermodynamic data. These values are derived either directly from the input data or
        calculated based on the provided configurations.

        Parameters
        ----------
        data : object
            The data object containing parameters, species, and reactions information
            required for initializing the kinetic model.

        Attributes
        ----------
        writer : Writer
            A writer instance used for logging messages.
        data : object
            A deep copy of the input data object.
        k_rate : None or float
            The kinetic rate, to be calculated or initialized later.
        electronic_part : None or float
            A placeholder for electronic contribution, to be set later.
        parameters : object
            The parameters from the input data object, used to configure
            properties of the kinetic model.
        species : object
            The species information from the input data object.
        reactions : object
            The reactions information from the input data object.
        operation : object
            A reference to the parameters for operations.
        electrode : float
            The electrode property, set to 1.0 for cathode or -1.0 for anode
            based on the input data.
        pre_exp : float
            The pre-exponential factor for rate constants, calculated based
            on the input data.
        experimental_part : numpy.ndarray
            A 2D array containing experimental rate constants for both
            forward and backward reactions, defaulted to ones if no experimental data is provided.
        thermochemical_part : numpy.ndarray
            A 2D array containing thermodynamic contributions to forward
            and backward reaction constants, initialized to zeros if no thermodynamic
            configurations are available.
        g_activation : numpy.ndarray or None
            Free energy of activation for reactions, calculated if thermochemical data is used.
        g_formation : numpy.ndarray or None
            Free energy of formation for adsorbed species, derived if thermochemical data is used.
        dg_reaction : numpy.ndarray or None
            Free energy of reaction, derived from thermodynamic contributions if configured.
        """
        self.data = copy.deepcopy(data)
        if writer is None:
            writer = Writer(
                log_file="melektrodica.log", log_directory=self.data.directory
            )
        writer.message(f"*** Kpynetic :  ***")

        super().__init__(self.data)
        self.k_rate = None
        self.electronic_part = None
        self.data = data
        self.parameters = data.parameters
        self.species = data.species
        self.reactions = data.reactions
        self.operation = data.parameters

        self.electrode = 1.0
        if not self.parameters.anode:
            self.electrode = -1.0

        # Pre-exponential
        self.pre_exp = self.parameters.pre_exponential
        if self.parameters.js:
            self.pre_exp = self.parameters.js_value / F
        if self.parameters.tst:
            self.pre_exp = (
                    self.parameters.kappa
                    * k_B
                    * self.parameters.temperature ** self.parameters.m
                    / h
            )

        # Experimental
        self.experimental_part = np.ones((2, len(self.reactions.list)))
        if self.parameters.experimental:
            self.experimental_part = RateConstants.experimental(
                self.reactions.k_f, self.reactions.k_b
            )

        # Thermodynamic
        self.thermochemical_part = np.zeros((2, len(self.reactions.list)))
        if self.parameters.thermochemical:
            self.g_activation = self.reactions.ga
            self.g_formation = self.species.g_formation_ads
            self.dg_reaction = FreeEnergy.reaction(
                self.reactions.upsilon_a, self.g_formation
            )
            self.thermochemical_part = RateConstants.thermochemical(
                self.g_activation, self.g_formation, self.reactions.upsilon_a
            )
            # if self.parameters.dg_reaction:
            #    self.dg_reaction = self.reactions.dg_reaction
            # elif self.parameters.g_formation:

    def foverpotential(
            self,
            potential: float,
            c_reactants: np.ndarray,
            c_products: np.ndarray,
            theta: np.ndarray,
    ) -> ndarray:
        """
        Calculate the overpotential and reaction rate within an electrochemical system.

        This function computes the overpotential and reaction rate for a given
        electrochemical reaction, based on the provided potential, concentrations
        of reactants and products, and surface coverage. Overpotential is calculated
        based on the electrode potential and reaction-specific parameters, including
        number of electrons transferred and symmetry factor. The reaction rate is
        determined using a combination of rate constants accounting for thermochemical,
        electronic, and experimental contributions.

        Parameters
        ----------
        potential : float
            The applied potential of the electrode.
        c_reactants : float
            The concentration of reactants in the system.
        c_products : float
            The concentration of products in the system.
        theta : float
            Surface coverage of the reaction intermediate.

        Notes
        -----
        The calculated overpotential `eta` is used in determining the
        electronic component of the reaction rate constant. This constant,
        along with the pre-exponential and experimental parts, is combined
        to compute the full rate constant. The rate constant is then used
        to determine the overall reaction rate `v`.
        """

        eta = potential
        self.electronic_part = self.electrode * RateConstants.electronic(
            eta, self.reactions.ne, self.reactions.beta
        )
        self.k_rate = self.constant(
            pre_exponential=self.pre_exp,
            experimental=self.experimental_part,
            thermochemical=self.thermochemical_part,
            electronic=self.electronic_part,
        )
        self.nu = self.rate(
            self.k_rate, c_reactants, c_products, theta, self.reactions.upsilon
        )

    def get_argument(self, potential: float):
        """
        Calculates and assigns the electronic part using the specified potential and reactions'
        parameters. Then computes and sets the necessary constant components required for
        further calculations. Finally, returns the calculated argument.

        Parameters
        ----------
        potential : float
            The potential value used to calculate the electronic part based on the
            electrode, number of electrons in the reaction, and the beta parameter.

        Returns
        -------
        float
            The computed argument based on the calculated electronic and other pre-defined
            constant components.
        """
        self.electronic_part = self.electrode * RateConstants.electronic(
            potential, self.reactions.ne, self.reactions.beta
        )
        self.constant(
            pre_exponential=1,
            experimental=1,
            thermochemical=self.thermochemical_part,
            electronic=self.electronic_part,
        )
        return self.argument

    def current(
            self,
            potential: float,
            c_reactants: np.ndarray,
            c_products: np.ndarray,
            theta: np.ndarray,
    ) -> ndarray:
        """
        Calculates the kinetic current density.

        .. math::

            J = F\\sum_i n_i \\nu_i

        This method models the relationship between the potential, concentrations of reactants and
        products, and other factors to compute the net current. The computation utilizes the
        overpotential and reaction mechanisms to derive the output.

        Parameters
        ----------
        potential : float
            The electrical potential in volts influencing the reaction.

        c_reactants : numpy.ndarray
            Concentrations of reactant species, represented as an array of values.

        c_products : numpy.ndarray
            Concentrations of product species, represented as an array of values.

        theta : float
            Surface coverage or activity parameter related to the reaction kinetics.

        Returns
        -------
        float
            The total electric current resulting from the electrochemical reactions.
        """
        self.foverpotential(potential, c_reactants, c_products, theta)
        return np.dot(self.reactions.ne, self.nu) * F

    def dcdt(self, rate, upsilon: np.ndarray) -> np.ndarray:
        """
        Computes the dot product of a given rate and upsilon values.

        .. math::
            \\frac{\\partial c_j}{\\partial t} =\\sum_i \\upsilon_{ij} \\nu_i

        .. math::
            \\frac{\\partial \\theta_j}{\\partial t} =\\sum_i \\upsilon_{ij} \\nu_i

        This method calculates the dot product of a rate value with an array
        of upsilon values using the numpy dot function. It is designed to
        return a numpy array representing the result of the computation.

        Parameters
        ----------
        rate : float
            A numerical value that represents the rate for the computation.
        upsilon : np.ndarray
            A numpy array containing the upsilon values to be used in the dot
            product computation.

        Returns
        -------
        np.ndarray
            A numpy array containing the result of the dot product computation.

        """
        return np.dot(rate, upsilon)

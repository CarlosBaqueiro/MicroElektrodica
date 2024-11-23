"""

    μElektrodica (Uxmal, version 1.0.0)
        A Python Tool for Modeling Microkinetic Electrocatalytic Reactions
        Copyright (C) 2024 C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        All rights reserved.

        Calculator class

"""

import numpy as np
from scipy.optimize import fsolve
from .Kpynetic import Kpynetic

# for debugging
#import sys
# sys.exit()
#from .Tools import showme

class BaseConcentration:
    """
    Represents a base class for concentration calculations in a chemical system.

    This class is designed to perform concentration and related state variable
    calculations for a given system defined by the parameters, species, and
    reactions.

    :ivar data: Data containing parameters, species, and reactions of the system.
    :ivar operation: Parameters related to the operation of the system.
    :ivar species: Species involved in the chemical reactions.
    :ivar reactions: Reactions occurring in the system.
    :ivar Kpy: Object containing specific methods for kinetic and potential calculations.
    :ivar c_reactants: Concentrations of reactants.
    :ivar c_products: Concentrations of products.
    :ivar theta: Coverages of adsorbed species.
    :ivar j: Array for current density calculations.
    :ivar fval: Array for function values from steady-state equations.
    """
    def __init__(self, data, Kpy):

        self.data = data
        self.operation = data.parameters
        self.species = data.species
        self.reactions = data.reactions
        self.Kpy = Kpy
        self.c_reactants = None
        self.c_products = None
        self.theta = None
        self.j = None
        self.fval = None



    def solver(self):
        """
        Solves the steady-state equations for a series of potentials and updates the concentrations of reactants,
        products, and adsorbed species, as well as the reaction current and function values.

        The method iterates over the potentials in the operation and uses a nonlinear solver `fsolve` to find the
        steady-state solution. It then updates the corresponding attributes with the calculated values.

        :return: Updated instance of the caller object with calculated values.
        :rtype: self
        """
        self.c_reactants = np.zeros((len(self.operation.potential), len(self.species.reactants)))
        self.c_products = np.zeros((len(self.operation.potential), len(self.species.products)))
        self.theta = np.zeros((len(self.operation.potential), len(self.species.adsorbed)))
        self.j = np.zeros(len(self.operation.potential))
        self.fval, initio = self.initialize()

        for i in range(len(self.operation.potential)):
            potential = self.operation.potential[i]
            solution = fsolve(self.steady_state, initio, args=potential, xtol=1e-12, maxfev=2000)
            self.c_reactants[i], self.c_products[i], self.theta[i] = self.unzip_variables(solution)
            self.fval[i] = self.steady_state(solution, potential)
            self.j[i] = self.current(solution, potential)
            initio = solution
        return self

    def initialize(self):
        """
        BaseClass serves as an abstract base class which requires subclasses
        to implement the `initialize` method. This class is designed to be
        inherited by other classes that provide specific implementations of
        initialization procedures.

        :raises NotImplementedError: when the `initialize` method is called
            on an instance of the BaseClass or any subclass that has not
            overridden this method.
        """
        raise NotImplementedError("The method initialize must be implemented by the subclass")

    def unzip_variables(self, variables):
        """
        Unzips a list of variables. Each subclass must implement this method
        to provide specific functionality for unzipping.

        :param variables: A list of variables to be unzipped.
        :type variables: list
        :raises NotImplementedError: If the method is not implemented by the subclass.
        :return: None
        """
        raise NotImplementedError("The method unzip_variables must be implemented by the subclass")

    def ride_hand_side(self, c_reactants, c_products, theta):
        """
        Compute the rate of change of reactants and products in a chemical reaction.

        This method calculates the right-hand side of a differential equation representing the
        rate of change of concentration of reactants and products in a chemical reaction
        system. The reaction can be characterized by specific parameters such as
        stoichiometric coefficients and kinetic parameters.

        :param c_reactants: Concentration of reactants.
        :type c_reactants: float
        :param c_products: Concentration of products.
        :type c_products: float
        :param theta: A parameter vector that may include rate constants and other
            kinetic parameters.
        :type theta: list[float]

        :return: Rate of change for reactants and products.
        :rtype: float
        """
        raise NotImplementedError("The method unzip_variables must be implemented by the subclass")

    def steady_state(self, variables, potential):
        """
        Calculates the steady-state concentration changes for reactants
        and products based on the given variables and potential.

        The method unwraps the variables into reactant concentrations,
        product concentrations, and an additional parameter theta. It then
        computes the right-hand side (rhs) of the equations using these
        concentrations. The potential function is updated with the
        concentrations and potential. Finally, it calculates the rate of
        change of concentrations and returns the difference between this
        rate and the rhs value.

        :param variables: List containing concentrations for reactants,
                          products and an additional parameter theta
                          in the order [c_reactants, c_products, theta].
        :type variables: list of float
        :param potential: A potential value affecting the concentration
                          changes.
        :type potential: float
        :return: Difference between the rate of change of concentrations
                 and the rhs value.
        :rtype: float
        """
        c_reactants, c_products, theta = self.unzip_variables(variables)
        rhs = self.ride_hand_side(c_reactants, c_products, theta)
        self.Kpy.potential_function(potential, c_reactants, c_products, theta)
        return self.Kpy.dcdt(self.Kpy.v, self.reactions.nux) - rhs

    def current(self, variables, potential):
        """
        Calculates the current based on the given variables and potential.

        This function utilizes the reaction-related concentrations and parameter
        theta from the extracted variables to compute the current by delegating
        the computation to the Kpy object.

        :param variables: A sequence containing concentrations of reactants,
            concentrations of products, and the parameter theta.
        :type variables: Sequence
        :param potential: The electric potential at which the current is calculated.
        :type potential: float
        :return: The calculated current.
        :rtype: float
        """
        c_reactants, c_products, theta = self.unzip_variables(variables)
        return self.Kpy.current(potential, c_reactants, c_products, theta)


class StaticConcentration(BaseConcentration):
    """
    Represents a static concentration model within a reaction system.

    This class initializes the concentrations of species within a reactive
    system and handles variables associated with the reaction process.
    It provides methods to initialize the system, unpack variables used in
    the reaction, and calculate the rate of change of the system.

    :ivar species: Contains the species involved in the reaction.
    :type species: Species
    :ivar operation: Holds operational parameters for the reaction.
    :type operation: Operation
    """
    def initialize(self):
        """
        Initializes the potential and adsorbed species arrays.

        :return: A tuple containing the initialized potential values array
        and the concatenated initial adsorbed species array.
        :rtype: tuple
        """
        fval = np.zeros((len(self.operation.potential), len(self.species.adsorbed)))
        theta0 = np.zeros(len(self.species.adsorbed))
        initio = np.concatenate([theta0])
        return fval, initio

    def unzip_variables(self, variables):
        """
        Unzips the provided variables and returns initial concentrations for reactants, products,
        and theta value from the species.

        :param variables: The variables to unzip.
        :type variables: any
        :return: A tuple containing the initial concentrations of reactants, initial concentrations of
                 products, and theta value.
        :rtype: tuple
        """
        c_reactants = self.species.c0_reactants
        c_products = self.species.c0_products
        theta = variables
        return c_reactants, c_products, theta

    def ride_hand_side(self, c_reactants, c_products, theta):
        """
        Compute the right-hand side of the reaction rate equation.

        This method calculates the right-hand side of the reaction rate equation given the
        concentrations of reactants and products, and a set of parameters.

        :param c_reactants: Concentrations of the reactants
        :type c_reactants: numpy.ndarray
        :param c_products: Concentrations of the products
        :type c_products: numpy.ndarray
        :param theta: Parameters of the reaction rate equation
        :type theta: numpy.ndarray
        :return: An array of zeros with the same length as the number of parameters in theta
        :rtype: numpy.ndarray
        """
        return np.zeros(len(theta))

class DynamicConcentration(BaseConcentration):
    """
    Handles dynamic concentration calculations for chemical species in a system.

    This class is responsible for initializing concentrations of reactants, products,
    and adsorbed species, and for unzipping and calculating the right-hand side of the
    differential equations governing the system.

    :ivar operation: Operation parameters containing information about potential, flow rate,
                     and cross-sectional area.
    :type operation: OperationParameters
    :ivar species: Species-related parameters including initial concentrations of reactants,
                   products, and adsorbed species.
    :type species: SpeciesParameters
    """
    def initialize(self):
        """
        Initializes the system by setting initial values for reactants, products, and adsorbed species.

        This function creates initial concentrations and initializes conditions for numerical
        simulations by setting up arrays of zeros and ones based on the lengths of reactants,
        products, and adsorbed species provided by the `species` attribute of the class instance.

        :return: A tuple containing the initial condition vectors.
        :rtype: tuple(np.ndarray, np.ndarray)
        """
        fval = np.zeros((len(self.operation.potential),
                         len(self.species.reactants) + len(self.species.products) + len(self.species.adsorbed)))
        c_reactants0 = np.ones(len(self.species.reactants))
        c_products0 = np.zeros(len(self.species.products))
        theta0 = np.zeros(len(self.species.adsorbed))
        initio = np.concatenate([c_reactants0, c_products0, theta0])
        return fval, initio

    def unzip_variables(self, variables):
        """
        Unzips a list of variables into reactants, products, and adsorbed species
        concentrations.

        :param variables: List of variables to unzip. The list is expected to be ordered such that
                          reactants are followed by products and then adsorbed species.
        :type variables: list
        :return: Tuple containing the concentrations of reactants, products, and adsorbed species.
        :rtype: tuple
        """
        c_reactants = variables[:len(self.species.reactants)]
        c_products = variables[len(self.species.reactants):-len(self.species.adsorbed)]
        theta = variables[-len(self.species.adsorbed):]
        return c_reactants, c_products, theta

    def ride_hand_side(self, c_reactants: np.ndarray, c_products: np.ndarray, theta: np.ndarray) -> np.ndarray:

        """
            Compute the right-hand side of the equation for the given reactant, product concentrations, and parameter theta.
            This function takes the concentrations of reactants and products and computes the right-hand side of the equation.
            It returns a concatenated numpy array containing the computed values for the reactants, products, and a zero-filled
            array of the same length as the theta parameter. This is typically used for simulations or mathematical modeling
            involving chemical species and their concentrations.

            :param c_reactants: Concentrations of reactants.
            :type c_reactants: np.ndarray
            :param c_products: Concentrations of products.
            :type c_products: np.ndarray
            :param theta: Parameter array theta, often related to system parameters or conditions.
            :type theta: np.ndarray

            :return: A concatenated array containing the computed values for reactants, products, and a zero-filled array
                     of the same length as theta.
            :rtype: np.ndarray
            """
        return np.concatenate([(c_reactants - self.species.c0_reactants) * self.operation.Fv / self.operation.Ac,
                               (c_products - self.species.c0_products) * self.operation.Fv / self.operation.Ac,
                               np.zeros(len(theta))])

class Calculator:
    """
    A class used to perform calculations based on provided data and chosen strategy.

    This class initializes with given data and selects a strategy for concentration
    calculation based on whether 'cstr' is set in the operation parameters. It then
    performs the calculations and stores the results.

    :ivar data: The initial data provided for calculations.
    :type data: dict
    :ivar operation: The operational parameters derived from the initial data.
    :type operation: dict
    :ivar Kpy: An instance of the Kpynetic class used for kinetic calculations.
    :type Kpy: Kpynetic
    :ivar strategy: The concentration calculation strategy chosen based on
        operational parameters.
    :type strategy: DynamicConcentration or StaticConcentration
    :ivar results: The results obtained from applying the chosen strategy's solver.
    :type results: dict
    """
    def __init__(self, data):
        """
        This class initializes with data and sets up a strategy for solving based
        on the given data parameters. It dynamically selects between dynamic and
        static concentration strategies, solving and storing the results.

        :param data: Input data containing parameters required for initialization.
        :type data: DataType

        :ivar data: Stores the input data.
        :ivar operation: Extracted parameters from the input data.
        :ivar Kpy: Instance of Kpynetic class initialized with data.
        :ivar strategy: Selected strategy based on the presence of cstr in data
                        parameters. Either DynamicConcentration or
                        StaticConcentration.
        :ivar results: Results obtained after solving using the selected strategy.
        """
        self.data = data
        self.operation = data.parameters

        self.Kpy = Kpynetic(self.data)
        if self.operation.cstr:
            self.strategy = DynamicConcentration(self.data, self.Kpy)
        else:
            self.strategy = StaticConcentration(self.data, self.Kpy)

        self.results = self.strategy.solver()



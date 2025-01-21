#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Calculator class

"""
import copy
import warnings
import numpy as np
from scipy.optimize import fsolve

from .kpynetic import Kpynetic
from .writer import Writer


# for debugging
# import sys
# sys.exit()
# from .Tools import showme


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

    def __init__(self, kpy):
        self.Kpy = kpy
        self.data = self.Kpy.data
        self.operation = self.data.parameters
        self.potential = self.operation.potential
        self.species = self.data.species
        self.reactions = self.data.reactions

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
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

        self.c_reactants = np.zeros(
            (len(self.operation.potential), len(self.species.reactants))
        )
        self.c_products = np.zeros(
            (len(self.operation.potential), len(self.species.products))
        )
        self.theta = np.zeros(
            (len(self.operation.potential), len(self.species.adsorbed))
        )
        self.j = np.zeros(len(self.potential))
        self.fval, initio = self.initialize()
        for i, potential in enumerate(self.operation.potential):
            solution = fsolve(
                self.steady_state, initio, args=potential, xtol=1e-12, maxfev=2000
            )
            self.c_reactants[i], self.c_products[i], self.theta[i] = (
                self.unzip_variables(solution)
            )
            self.fval[i] = self.steady_state(solution, potential)
            self.j[i] = self.current(solution, potential)
            initio = solution

        for warning in w:
            if issubclass(warning.category, RuntimeWarning):
                self.Kpy.writer.logger.error(f"WARNING at potential {potential}: {warning.message}")
                if "The iteration is not making good progress" in str(warning.message):
                    raise RuntimeError(f"Convergence failed at potential {potential}")
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
        raise NotImplementedError(
            "The method initialize must be implemented by the subclass"
        )

    def unzip_variables(self, variables):
        """
        Unzips a list of variables. Each subclass must implement this method
        to provide specific functionality for unzipping.

        :param variables: A list of variables to be unzipped.
        :type variables: list
        :raises NotImplementedError: If the method is not implemented by the subclass.
        :return: None
        """
        raise NotImplementedError(
            "The method unzip_variables must be implemented by the subclass"
        )

    def right_hand_side(self, c_reactants, c_products, theta):
        """
        Compute the right-hand side of the reaction rate equation.

        This method calculates the right-hand side of the reaction rate equation given the
        concentrations of reactants and products, and a set of parameters.

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
        raise NotImplementedError(
            "The method unzip_variables must be implemented by the subclass"
        )

    def steady_state(self, variables, potential):
        """
        Calculates the steady-state concentration changes for reactants
        and products based on the given variables and potential.

        The method unwraps the variables into reactant concentrations,
        product concentrations, and an additional parameter theta. It then
        computes the right-hand side (rhs) of the equations using these
        concentrations. The potential function is updated with the
        concentrations and potential. Finally, it calculates the rate
        change of concentrations and returns the difference between this
        rate and the rhs value.

        :param variables: List containing concentrations for reactants,
                          products and an additional parameter theta
                          in the order [c_reactants, c_products, theta].
        :type variables: list of float
        :param potential: A potential value affecting the concentration
                          changes.
        :type potential: float
        :return: Difference between the rate change of concentrations
                 and the rhs value.
        :rtype: float
        """
        c_reactants, c_products, theta = self.unzip_variables(variables)
        rhs = self.right_hand_side(c_reactants, c_products, theta)
        self.Kpy.foverpotential(potential, c_reactants, c_products, theta)
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
    the reaction, and calculate the rate change of the system.

    :ivar species: Contains the species involved in the reaction.
    :type species: Species
    :ivar operation: Holds operational parameters for the reaction.
    :type operation: Operation
    """

    def initialize(self):
        """
        Initializes the system by creating initial values and arrays for the simulation.

        This function prepares initial values of the state variables required for a
        simulation. It initializes two arrays: one for storing function values and
        another with initial condition values derived from the given data structure.
        The initial states are related to adsorption species and potential operations.

        Returns
        -------
        tuple
            A tuple containing:
            - fval : ndarray
                An array initialized with zeros, with shape based on the lengths
                of potential operations and adsorbed species.
            - initio : ndarray
                A 1-dimensional array created by concatenating initial values
                pertaining to adsorbed species.
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

    def right_hand_side(self, c_reactants, c_products, theta):
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

        fval = np.zeros(
            (
                len(self.operation.potential),
                len(self.species.reactants)
                + len(self.species.products)
                + len(self.species.adsorbed),
            )
        )
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
        c_reactants = variables[: len(self.species.reactants)]
        c_products = variables[
                     len(self.species.reactants): -len(self.species.adsorbed)
                     ]
        theta = variables[-len(self.species.adsorbed):]
        return c_reactants, c_products, theta

    def right_hand_side(
            self, c_reactants: np.ndarray, c_products: np.ndarray, theta: np.ndarray
    ) -> np.ndarray:
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
        return np.concatenate(
            [
                (c_reactants - self.species.c0_reactants)
                * self.operation.Fv
                / self.operation.Ac,
                (c_products - self.species.c0_products)
                * self.operation.Fv
                / self.operation.Ac,
                np.zeros(len(theta)),
            ]
        )


class Calculator:
    """
    Represents a computational tool for executing operations, simulations, and
    evaluations within a specified system.

    Provides mechanisms for handling system data, operational parameters,
    and computational strategies, alongside logging features for tracking
    execution and results. Assigns specific calculation strategies dynamically
    based on the characteristics of the operational configuration.

    Attributes:
        name: str
            The name assigned to the calculator instance. If not provided, defaults to 'melek'.
        writer: Writer
            An instance of the Writer class used for logging and messaging.
        Kpy: <type>
            A deep copy of the input kpy object used for operations and computations.
        data: <type>
            Holds the data contained within the Kpy object, representing input information for the calculations.
        operation: <type>
            References the `parameters` attribute of the data object and represents operational parameters.
        potential: <type>
            References the `potential` attribute within the operational parameters.
        species: <type>
            Represents all species involved in the system, extracted from the data object.
        reactions: <type>
            Represents all reactions within the system, extracted from the data object.
        strategy: DynamicConcentration or StaticConcentration
            Chooses between dynamic or static concentration strategies based on the cstr attribute in operation.
        results: <type>
            Stores the output of the solver executed from the selected strategy.

    Raises:
        ValueError
            Raised if the results from the strategy solver contain negative values in `theta`.
    """

    def __init__(self, kpy, name=None):
        """
        Initializes the Calculator class and sets up default attributes and related objects required for computational processes.

        Attributes:
            name: str
                The name assigned to the calculator instance. If not provided, defaults to 'melek'.
            writer: Writer
                An instance of the Writer class used for logging and messaging.
            Kpy: <type>
                A deep copy of the input kpy object used for operations and computations.
            data: <type>
                Holds the data contained within the Kpy object, representing input information for the calculations.
            operation: <type>
                References the `parameters` attribute of the data object and represents operational parameters.
            potential: <type>
                References the `potential` attribute within the operational parameters.
            species: <type>
                Represents all species involved in the system, extracted from the data object.
            reactions: <type>
                Represents all reactions within the system, extracted from the data object.
            strategy: DynamicConcentration or StaticConcentration
                Chooses between dynamic or static concentration strategies based on the cstr attribute in operation.
            results: <type>
                Stores the output of the solver executed from the selected strategy.

        Raises:
            ValueError
                Raised if the results from the strategy solver contain negative values in `theta`.
        """
        if name is None:
            self.name = 'melek'
        else:
            self.name = name

        self.writer = Writer()
        self.writer.message(f"*** Calculator : {self.name}  ***")

        self.Kpy = copy.deepcopy(kpy)
        self.data = self.Kpy.data
        self.operation = self.data.parameters
        self.potential = self.operation.potential
        self.species = self.data.species
        self.reactions = self.data.reactions

        if self.operation.cstr:
            self.strategy = DynamicConcentration(self.Kpy)
        else:
            self.strategy = StaticConcentration(self.Kpy)

        #def strategy_solver(self):
        self.results = self.strategy.solver()
        if np.any(self.results.theta < 0):
            self.writer.logger.error("Solution contains negative values")

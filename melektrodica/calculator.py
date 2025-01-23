#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    μElektrodica © 2025
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
    Manages the initialization and storage of chemical reaction model data.

    This class is designed to initialize and organize key parameters and
    data for simulation or analysis of chemical reaction models. It sets
    up critical attributes such as operation parameters, potential, species,
    and reactions, while also providing placeholders for computed values.

    Attributes
    ----------
    Kpy : object
        The reference to the input `kpy` object provided during initialization.
    data : object
        The data attribute of `kpy`, containing necessary information for
        reaction modeling.
    operation : object
        A sub-attribute of `data`, representing the operational parameters.
    potential : object
        Extracted from `operation`, represents the chemical or model potential.
    species : object
        Extracted from `data`, contains the chemical species present in the model.
    reactions : object
        Extracted from `data`, representing the set of reactions in the chemical
        model.
    c_reactants : None or other
        Placeholder for computed reactant concentrations.
    c_products : None or other
        Placeholder for computed product concentrations.
    theta : None or other
        Placeholder for a computed parameter (e.g., coverage fraction).
    j : None or other
        Placeholder for computed current density or rate parameter.
    fval : None or other
        Placeholder for a computed function value during calculations
        (e.g., objective function value).
    """

    def __init__(self, kpy):
        """
        Manages the initialization and storage of chemical reaction model data.

        This class is designed to initialize and organize key parameters and data
        for simulation or analysis of chemical reaction models. It sets up critical
        attributes such as operation parameters, potential, species, and reactions,
        while also providing placeholders for computed values.

        Parameters
        ----------
        kpy : object
            The input object from which the chemical reaction data and parameters
            are extracted. This must contain `data` with `parameters`, `potential`,
            `species`, and `reactions` attributes.

        Attributes
        ----------
        Kpy : object
            The reference to the input `kpy` object provided during initialization.
        data : object
            The data attribute of `kpy`, containing necessary information for reaction
            modeling.
        operation : object
            A sub-attribute of `data`, representing the operational parameters.
        potential : object
            Extracted from `operation`, represents the chemical or model potential.
        species : object
            Extracted from `data`, contains the chemical species present in the model.
        reactions : object
            Extracted from `data`, representing the set of reactions in the chemical
            model.
        c_reactants : None or other
            Placeholder for computed reactant concentrations.
        c_products : None or other
            Placeholder for computed product concentrations.
        theta : None or other
            Placeholder for a computed parameter (e.g., coverage fraction).
        j : None or other
            Placeholder for computed current density or rate parameter.
        fval : None or other
            Placeholder for a computed function value during calculations (e.g.,
            objective function value).
        """
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
            solver(self)

            Solves a system of equations for steady-state reaction kinetics and computes
            reactant, product, and adsorbed species concentrations as well as the
            current density for a range of applied potentials. The method utilizes
            numerical root-finding techniques to achieve steady-state solutions.

            During the solving process, runtime warnings are captured and logged for
            debugging purposes. If convergence issues occur, an exception is raised
            for the specific potential value.

            Returns
            -------
            self : object
                The instance of the class with updated attributes for steady-state
                reactant, product, adsorbed species concentrations, computed current
                densities, and other intermediate results.
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
                self.steady_state, initio, args=potential, xtol=1e-9, maxfev=2000
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
        Raises
        ------
        NotImplementedError
            This exception is raised when the method `initialize` is not implemented
            by a subclass. Subclasses inheriting this method must override and
            implement their own logic.
        """
        raise NotImplementedError(
            "The method initialize must be implemented by the subclass"
        )

    def unzip_variables(self, variables):
        """
        Unzips a collection of variable pairs into two separate lists.

        This method is intended to be overridden by subclasses to provide the specific
        implementation for unzipping the given variables. It takes a collection of
        paired variables and separates them into two distinct lists or structures.

        Parameters
        ----------
        variables : list of tuple
            A collection of paired variables, where each pair is represented as a tuple.
            The method expects the input to be iterable.

        Raises
        ------
        NotImplementedError
            If this method is not overridden by a subclass and is called directly,
            this exception is raised to indicate that the implementation is missing.
        """
        raise NotImplementedError(
            "The method unzip_variables must be implemented by the subclass"
        )

    def right_hand_side(self, c_reactants, c_products, theta):
        """
        Computes the right-hand side of a system of ordinary differential equations (ODEs).

        This method calculates the derivative of concentrations with respect to time, using
        the provided reactant concentrations, product concentrations, and kinetic parameters.
        The actual implementation must be provided in a subclass by overriding this method.

        Parameters
        ----------
        c_reactants : Any
            The concentrations of the reactants in the system. The specific data type and
            structure depend on the implementation and the problem being modeled.
        c_products : Any
            The concentrations of the products in the system. The specific data type and
            structure depend on the implementation and the problem being modeled.
        theta : Any
            The vector or parameters representing the kinetic coefficients or other
            parameters necessary for computing the reaction rates.

        Raises
        ------
        NotImplementedError
            If the method is not overridden in a subclass.
        """
        raise NotImplementedError(
            "The method unzip_variables must be implemented by the subclass"
        )

    def steady_state(self, variables, potential):
        """
        Computes the steady-state properties of a reaction system given the input variables and
        potential. This function evaluates the reaction dynamics by decomposing the input
        variables, computing the right-hand side (RHS) of the reaction system, evaluating the
        overpotential, and determining the change in concentrations and other state properties.

        Parameters
        ----------
        variables : Any
            A set of state variables for the reaction system, which include the concentrations
            of reactants, products, and the state variable `theta`.

        potential : Any
            The potential applied to the reaction system which is utilized for calculating the
            overpotential.

        Returns
        -------
        Any
            The computed rate of change (`dcdt`) for the concentrations and other reaction
            variables, determined by subtracting the right-hand side from the calculated rate
            expressions.
        """

        c_reactants, c_products, theta = self.unzip_variables(variables)
        rhs = self.right_hand_side(c_reactants, c_products, theta)
        self.Kpy.foverpotential(potential, c_reactants, c_products, theta)
        return self.Kpy.dcdt(self.Kpy.v, self.reactions.upsilonx) - rhs

    def current(self, variables, potential):
        """
        Calculates the current for a given set of variables and potential using
        an underlying kinetic model.

        This function extracts the necessary parameters from the input variables
        and processes the potential and concentrations of reactants and products.
        The current is then computed using a predefined kinetic model instance.

        Parameters
        ----------
        variables : Any
            The variables containing information about reactant concentrations,
            product concentrations, and additional parameters required for
            current calculation.
        potential : Any
            The potential at which the current is to be calculated. This is a key
            input for the kinetic model.

        Returns
        -------
        Any
            The calculated current value based on the provided potential and
            extracted variables.
        """
        c_reactants, c_products, theta = self.unzip_variables(variables)
        return self.Kpy.current(potential, c_reactants, c_products, theta)


class StaticConcentration(BaseConcentration):
    """
    Simulates the static concentration within a computational model.

    This class is designed to handle the initialization, transformation,
    and computation of parameters and variables in modeling adsorption
    species and their interactions with potentials. It includes methods
    to set up the initial state of the system, extract variables into
    reactants and products, and compute specific reaction rate equations.

    Attributes
    ----------
    species : Species
        Contains details related to adsorbed species, including initial
        concentrations.
    operation : Operation
        Represents the specific potentials and operational dimensions
        involved in the computation process.
    """

    def initialize(self):
        """
        initialize(self)

        Initializes and sets up the necessary initial state for a simulation or computation
        process involving adsorption species and potential operations. Specifically, this
        method creates zero matrices and vectors based on the dimensions of the adsorption
        species and the potentials involved, preparing for the subsequent computations.

        Returns
        -------
        tuple of (ndarray, ndarray)
            A tuple containing:
            - fval : ndarray
                A 2D array of zeros with shape corresponding to the number of potentials
                and the number of adsorbed species.
            - initio : ndarray
                A 1D array of zeros concatenated from adsorption species initial states.
        """

        fval = np.zeros((len(self.operation.potential), len(self.species.adsorbed)))
        theta0 = np.zeros(len(self.species.adsorbed))
        initio = np.concatenate([theta0])
        return fval, initio

    def unzip_variables(self, variables):
        """
        Extracts the initial concentrations for reactants and products along with a variable, theta.

        This function unpacks a set of variables into its corresponding initial concentrations
        for reactants and products, and a specific variable referred to as theta. These
        concentrations are determined using the associated species attributes of the object.

        Parameters
        ----------
        variables : Any
            A variable or set of variables that contains the necessary components
            to unpack into `theta`.

        Returns
        -------
        tuple
            A tuple containing:
            - c_reactants: The initial concentrations for reactants.
            - c_products: The initial concentrations for products.
            - theta: The unpacked variable from the input variables.

        """

        c_reactants = self.species.c0_reactants
        c_products = self.species.c0_products
        theta = variables
        return c_reactants, c_products, theta

    def right_hand_side(self, c_reactants, c_products, theta):
        """
        Computes the right-hand side of a system of equations describing reaction dynamics.

        This method evaluates the rate equations for chemical reactions, which describe
        changes in concentrations of reactants and products over time based on reaction
        constants and reaction conditions. It returns the computed values of the system's
        right-hand side at a given state.

        Parameters
        ----------
        c_reactants : np.ndarray
            Array of concentrations of reactant species.
        c_products : np.ndarray
            Array of concentrations of product species.
        theta : np.ndarray
            Array of parameters (e.g., rate constants) for the reactions.

        Returns
        -------
        np.ndarray
            Array of evaluated right-hand side values corresponding to each reaction.
        """
        return np.zeros(len(theta))


class DynamicConcentration(BaseConcentration):
    """
    Handles dynamic concentration calculations for chemical species during simulations.

    This class is designed to model and compute the dynamic concentrations of reactants,
    products, and adsorbed species. It provides methods for initializing simulation conditions,
    unpacking concentration arrays, and computing the right-hand side of the governing
    differential equations for the modeled system.

    Attributes
    ----------
    species : object
        Object containing information about reactants, products, and adsorbed species. It
        should provide attributes like `reactants`, `products`, `adsorbed`, `c0_reactants`,
        and `c0_products`.
    operation : object
        Object containing operational parameters such as potential, fluid velocity (`Fv`),
        and cross-sectional area (`Ac`), which are used for system calculations.
    """

    def initialize(self):
        """
        Initializes and returns the state variables and a matrix of zeros for further computational
        operations within a chemical kinetic simulation. This function primarily prepares
        initial concentrations for reactants, products, and adsorbed species, and a zero-initialized
        matrix associated with potential data and species being modeled.

        Returns
        -------
        fval : numpy.ndarray
            A zero-initialized matrix with dimensions corresponding to the number of
            potential values in the `operation.potential` attribute across the
            summed count of reactant, product, and adsorbed species in the `species`
            attribute.

        initio : numpy.ndarray
            An array containing the concatenated initial concentration values for
            reactants (set to ones), products (set to zeros), and adsorbed species
            (set to zeros).
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
        Unzips a list of variables into separate components representing reactants,
        products, and adsorbed species concentrations.

        The method splits the input `variables` array into three distinct components:
        concentrations of reactants (`c_reactants`), concentrations of products
        (`c_products`), and the fractional surface coverage of adsorbed species
        (`theta`). The split is determined by the number of species categorized into
        reactants, products, and adsorbed species within the `self.species` object.

        Parameters
        ----------
        variables : list
            A one-dimensional list of variables representing the concentrations
            of reactants, products, and the surface coverage of adsorbed species.

        Returns
        -------
        tuple
            A tuple containing three elements:
            - c_reactants : list
                Concentrations of reactant species.
            - c_products : list
                Concentrations of product species.
            - theta : list
                Fractional surface coverage of adsorbed species.
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
        Computes the right-hand side of the set of differential equations governing the
        reaction system. This involves combining contributions from the reactants,
        products, and additional state variables (`theta`), considering the changes
        relative to their initial concentrations and system operational parameters.

        Parameters
        ----------
        c_reactants : np.ndarray
            Concentrations of the reactant species in the system.
        c_products : np.ndarray
            Concentrations of the product species in the system.
        theta : np.ndarray
            Additional state variables representing system-specific parameters or
            conditions.

        Returns
        -------
        np.ndarray
            A combined array representing the contributions of reactants, products, and
            additional state variables to the right-hand side of the equations.
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
    Represents a computational framework for processing data, strategies and calculations
    based on the input system parameters and settings.

    Attributes
    ----------
    name : str
        The name assigned to the calculator instance. Defaults to 'melek' if not specified.
    writer : Writer
        An instance of the Writer class used for logging and messaging functionalities.
    Kpy : deepcopy of kpy
        A deep copy of the input kpy object used for operations and computations.
    data : type inferred from kpy.data
        Contains the input information extracted from the Kpy object for calculations.
    operation : type inferred from data.parameters
        Represents the operational parameters extracted from the data object.
    potential : type inferred from operation.potential
        Refers to the `potential` parameter within the operational attributes.
    species : type inferred from data.species
        Denotes all species involved in the system, retrieved from the data object.
    reactions : type inferred from data.reactions
        Denotes all reactions in the current system, retrieved from the data object.
    strategy : DynamicConcentration or StaticConcentration
        The computational strategy applied, either dynamic or static concentration,
        depending on the `cstr` setting in operation.
    results : type inferred from strategy.solver()
        The output generated by executing the solver of the defined strategy.

    Raises
    ------
    ValueError
        Raised if the results computed from the strategy solver contain negative values in `theta`.
    """

    def __init__(self, kpy, name=None):
        """
        Initializes a Calculator instance and sets it up to calculate based on the provided
        kpy data structure. Determines the operation type (dynamic or static concentration) and
        computes the solution for the system.

        Parameters
        ----------
        kpy : object
            The input data structure that contains all necessary information for constructing
            the Calculator, including data, parameters, species, and reactions.

        name : str, optional
            The name of the calculator. If not provided, defaults to 'melek'.

        Attributes
        ----------
        name : str
            The name of the calculator instance.

        writer : Writer
            Handles logging and messaging for the instance.

        Kpy : object
            A deepcopy of the provided `kpy` input to prevent modifications to the original
            object.

        data : object
            Extracted data object from `kpy` containing various system configurations.

        operation : object
            Parameters for the operation derived from the `data` object.

        potential : narray
            Potential derived from the operation's parameters.

        species : object
            Species list derived from the `data`.

        reactions : object
            Reaction information derived from the `data`.

        strategy : object
            Dynamic or static concentration calculation strategy determined based on whether
            CSTR operation is specified in the `operation`.

        results : object
            Solution results obtained by applying the selected strategy's solver method.
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

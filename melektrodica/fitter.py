"""

    μElektrodica© 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Fitter class

"""

import copy
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
from scipy.optimize import Bounds, differential_evolution
from .calculator import Calculator
from .writer import Writer


# for debugging
# import sys
# sys.exit()
# from .Tools import showme


class Fitter(Calculator):
    """
    A specialized class for energy fitting and optimization.

    This class extends the Calculator class and is designed to perform energy
    fitting tasks using given thermodynamic and kinetic data. It includes functionalities
    for optimizing reaction and formation energies based on experimental data and
    implements utilities such as error tracking and optimization bounds.

    The least-squares parameter estimation problem for determining kinetic parameters is formulated as:

    .. math::

            f(\\mathbf{p}) = \\min_{\\mathbf{p}} \\;\\frac{1}{2}
                                \\sum \left[ \frac{J_{\\text{exp}}(\\eta)
                                -
                                J_{\\text{model}}(\\eta,\\mathbf{p})}{J_{\\text{exp}}(\\eta)} \\right]^2

    Attributes
    ----------
    name : str
        Name identifier for the Fitter instance.
    writer : Writer
        A writer instance to log messages or outputs.
    Kpy : deep copy of the passed object
        A deep copy of the input Kpy object.
    data : data structure
        Input data from Kpy, holding reactions, species, and their properties.
    potential_data : deep copy of the passed object
        Copy of potential data for use in the fitting process.
    j_data : deep copy of the passed object
        Experimental current data used as a reference for fitting comparison.
    error_evolution : list
        List to store the evolution of the objective function value during optimization.
    bounds : Bounds
        Optimization bounds for reaction and formation energy variables.
    g_fit : optimization result
        Results of the energy fitting optimization process.
    ga_fit : ndarray
        Optimized reaction energies derived from the fitted results.
    gf_fit : ndarray
        Optimized formation energies derived from the fitted results.
    """

    def __init__(self, kpy, potential_data, j_data, name=None):
        """
        Initializes the Fitter object and sets up necessary attributes and optimization bounds
        using provided data. This class manages data for a fitting process and performs
        initialization for scaling factors and parameter bounds.

        Attributes
        ----------
        name : str
            The name identifier for the fitter instance. Default is 'melek' if not provided.
        writer : Writer
            Object used for logging messages during initialization.
        Kpy : object
            Deep copy of the provided `kpy` object containing data for the fitter.
        data : object
            Data object extracted from `Kpy`, containing relevant fitting information.
        potential_data : object
            Deep copy of the provided potential data, which is used to define the system's
            parameters.
        j_data : object
            Deep copy of the provided journal data, which is potentially used for calibrations.
        error_evolution : list
            A list that will store the evolution of error during fitting. Empty at initialization.
        bounds : Bounds
            Optimization bounds specified for the fitting process, calculated based on lower
            and upper limits derived from the system's reactions and species formation energies.
        g_fit : OptimizeResult
            Optimization result containing the best-fit energy parameter values as a result
            of the fitting process.
        ga_fit : ndarray
            Fitted energy parameters related to reactions.
        gf_fit : ndarray
            Fitted energy parameters related to species formation energies.

        Parameters
        ----------
        kpy : object
            Input object containing data necessary for initializing the fitter, which includes
            system parameters and reactions.
        potential_data : object
            Potential data used for setting up the system's parameters in the fitting process.
        j_data : object
            Journal data used for additional calibration in the fitting process, if required.
        name : str, optional
            An optional name identifier for the fitter. If not provided, the default name 'melek'
            will be used.

        """
        if name is None:
            self.name = "melek"
        else:
            self.name = name

        self.writer = Writer()
        self.writer.message(f"*** Fitter : {self.name}  ***")

        self.Kpy = copy.deepcopy(kpy)
        self.data = self.Kpy.data
        self.potential_data = copy.deepcopy(potential_data)
        self.j_data = copy.deepcopy(j_data)
        self.data.parameters.potential = self.potential_data
        super().__init__(self.Kpy)
        self.error_evolution = []

        g0 = np.concatenate([self.data.reactions.ga, self.data.species.g_formation_ads])
        g_lb = g0 * 0.8  # Lower limit (80% of the initial value)
        g_ub = g0 * 1.2  # Upper limit (120% of the initial value)
        g_lb, g_ub = np.minimum(g_lb, g_ub), np.maximum(g_lb, g_ub)
        self.bounds = Bounds(lb=g_lb, ub=g_ub, keep_feasible=True)
        self.g_fit = self.fit_energies()
        self.ga_fit, self.gf_fit = self.unziper(self.g_fit.x)

    def fit_energies(self):
        """
        Optimize a given objective function using the Differential Evolution algorithm.

        This method utilizes the `differential_evolution` optimization algorithm to fit
        the provided objective function over specified bounds. The optimization process
        is configured using a variety of settings such as mutation, recombination, and
        population size, among others. If an error occurs during the optimization, it
        is caught and logged, and the method safely returns `None`.

        Returns
        -------
        OptimizeResult or None
            The result of the differential evolution optimization process if
            successful; otherwise, returns `None`.

        Raises
        ------
        Exception
            Captures any unexpected errors during the optimization process, logs the
            error message, and safely returns `None`.

        Notes
        -----
        This method is tailored for scenarios that involve energy-based optimization
        tasks. It includes a callback function to handle and display error evolution
        during the optimization process. The Differential Evolution algorithm explores
        the search space via a population-based stochastic approach, which is effective
        for nonlinear and non-differentiable functions.
        """
        try:
            g_fit = differential_evolution(
                func=self.object,
                bounds=self.bounds,
                strategy="best1bin",
                popsize=15,
                tol=1e-3,
                mutation=(0.5, 1.0),
                recombination=0.7,
                seed=None,
                disp=False,
                polish=True,
                init="latinhypercube",
                updating="immediate",
                callback=self.display_error_evolution,
            )
            return g_fit

        except Exception as e:
            print(f"Error during optimization: {e}")
            return None

    def object(self, *energies):
        """
        Calculate the fitting error between experimental and calculated currents.

        This function computes the squared error between experimentally obtained
        currents (`j_data`) and currents computed from a model (`j_fit`) based on the
        given `energies`. It tracks the error evolution over successive iterations,
        providing a measure of goodness of fit.

        Parameters
        ----------
        energies : tuple
            Energies used to calculate the modeled current (`j_fit`) using the
            `current_energies` method.

        Returns
        -------
        float
            The calculated fitting error. Returns infinity (`np.inf`) in cases where
            errors occur during calculation or zero values are present in
            experimental data (`j_data`).
        """
        try:
            j_fit = np.abs(self.current_energies(*energies))
            if np.any(self.j_data == 0):
                print("Warning: Encountered zero in calculated currents.")
                return np.inf

            # Calculate the squared error with respect to the experimental data
            error = np.sum(np.pow(self.j_data - j_fit, 2) / self.j_data)
            # print(f"Error: {error}")
            self.error_evolution.append(error)
            return error
        except Exception as e:
            print(f"Error in fitting calculation: {e}")
            return np.inf

    def unziper(self, variables):
        """
        Unzips the input tuple into two distinct parts: the first part corresponds
        to a subset aligned with the length of a specific attribute in the object,
        and the second part represents the remaining elements.

        The method primarily operates on tuples, extracting and splitting data
        based on pre-determined boundaries determined by the length of a specific
        attribute.

        Parameters
        ----------
        variables : tuple or any
            The input data to be unzipped. If it's a tuple, it is divided into two
            components based on the length of `self.data.reactions.list`. If not
            a tuple, it is processed as a single component.

        Returns
        -------
        tuple
            A tuple containing two elements:
            - The first element is aligned in size with the length of
              `self.data.reactions.list`.
            - The second element contains the rest of the items beyond that length.
        """
        if isinstance(variables, tuple):
            variables = variables[0]
        a = variables[: len(self.data.reactions.list)]
        f = variables[len(self.data.reactions.list) :]
        return a, f

    def current_energies(self, *energies):
        """
        Computes and updates the current energies provided as input, processes them
        with thermodynamic calculations, and computes results through a defined
        strategy solver.

        This method takes in variable energies, processes them by splitting into
        specific thermodynamic components, updates the related thermochemical
        attributes, and computes final results using a solver from the specified
        strategy. Errors are handled and printed if attributes are missing or
        unexpected issues occur during execution.

        Parameters
        ----------
        energies : tuple
            Variable length tuple representing the input energies that are
            processed and used for calculations within the method.

        Returns
        -------
        float
            The resultant value `j` from the calculations performed by the
            strategy solver after the energies are processed.

        Raises
        ------
        AttributeError
            Raised when an attribute accessed within the method is not found or
            defined.
        Exception
            Raised for any unexpected errors encountered during execution, with
            details of the error printed for debugging.
        """
        try:
            ga, gf = self.unziper(energies)
            # Update Kpy with thermodynamic values
            self.Kpy.thermochemical_part = self.Kpy.thermochemical(
                ga, gf, self.results.data.reactions.upsilon_a
            )
            self.results = self.strategy.solver()
            return self.results.j

        except AttributeError as e:
            print(f"Attribute error in current_energies: {e}")
            raise

        except Exception as e:
            print(f"Unexpected error in current_energies: {e}")
            raise

    def display_error_evolution(self, xk, convergence=0):
        """
        Displays the error evolution graph during the optimization process.

        This method creates a real-time plot that represents the evolution
        of the objective function values across iterations. It uses a logarithmic
        scale for the y-axis to better visualize the changes in error values and
        updates the plot dynamically as the optimization progresses. It is useful
        for monitoring the runtime optimization process and understanding the
        convergence behavior.

        Parameters
        ----------
        xk : Any
            Current parameter values during the optimization process. The parameter
            values are not used directly in this function but are part of the optimizer
            interface.
        convergence : int, optional
            Convergence tolerance or metric used in the optimization process. This
            parameter is not used directly in this plotting method. Default is 0.

        Notes
        -----
        This function uses `matplotlib.pyplot` for plotting and dynamically updates
        the graph using `plt.pause(0.01)` to provide a real-time visualization. To
        avoid overlapping outputs, the previous plot is cleared using
        `clear_output(wait=True)`.
        """
        clear_output(wait=True)
        plt.figure()
        plt.plot(self.error_evolution, label="Objective function value")
        plt.yscale("log")
        plt.xlabel("Iterations")
        plt.ylabel("Objective function value")
        plt.title("Runtime optimization")
        plt.legend()
        plt.show()
        plt.pause(0.01)

    def new_results(self, new_potential):
        """
        Updates the calculator with new potential values and calculates the respective
        thermochemical parameters utilizing the current state of the system.

        The method modifies the potential attribute in the parameters and recalculates
        the thermochemical part within the `Kpy`. After applying the thermodynamic
        calculations, it returns a new `Calculator` object initialized with the updated
        `Kpy` instance and the mode specified.

        Parameters
        ----------
        new_potential : float
            The new potential value to be assigned for recalculations.

        Returns
        -------
        Calculator
            A new Calculator instance initialized with the updated `Kpy` object
            and the mode set as 'Fitter'.
        """
        # Update Kpy with thermodynamic and potential values
        self.data.parameters.potential = new_potential
        self.Kpy.thermochemical_part = self.Kpy.thermochemical(
            self.ga_fit, self.gf_fit, self.results.data.reactions.upsilon_a
        )
        return Calculator(self.Kpy, 'Fitter')
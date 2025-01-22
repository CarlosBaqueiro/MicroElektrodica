"""

    μElektrodica © 2024
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
#import sys
# sys.exit()
# from .Tools import showme


class Fitter(Calculator):
    """
    Fitter class for optimizing reaction and species energy data to match experimental current values.
    Inherits from Calculator and integrates optimization methods.
    """

    def __init__(self, kpy, potential_data, j_data, name=None):
        if name is None:
            self.name = 'melek'
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
        Executes the optimization routine to fit energy data.
        Returns optimization results.
        """
        try:
            g_fit = differential_evolution(
                func=self.object,
                bounds=self.bounds,
                strategy='best1bin',
                popsize=15,
                tol=1e-3,
                mutation=(0.5, 1.0),
                recombination=0.7,
                seed=None,
                disp=False,
                polish=True,
                init='latinhypercube',
                updating='immediate',
                callback=self.display_error_evolution
            )
            return g_fit

        except Exception as e:
            print(f"Error during optimization: {e}")
            return None

    def object(self, *energies):
        """
        Computes the error between calculated and experimental values for given energies.
        """
        try:
            j_fit = np.abs(self.current_energies(*energies))
            if np.any(self.j_data == 0):
                print("Warning: Encountered zero in calculated currents.")
                return np.inf

            # Calculate the squared error with respect to the experimental data
            error = np.sum(np.pow(self.j_data - j_fit, 2) / self.j_data)
            #print(f"Error: {error}")
            self.error_evolution.append(error)
            return error
        except Exception as e:
            print(f"Error in fitting calculation: {e}")
            return np.inf

    def unziper(self, variables):
        """
        Splits a combined array of variables into reaction energies and formation energies.
        """
        if isinstance(variables, tuple):
            variables = variables[0]
        a = variables[:len(self.data.reactions.list)]
        f = variables[len(self.data.reactions.list):]
        return a, f

    def current_energies(self, *energies):
        """
        Calculates current energies for a given set of energies.
        """
        try:
            ga, gf = self.unziper(energies)
            # Update Kpy with thermodynamic values
            self.Kpy.thermodynamic_part = self.Kpy.thermodynamic(ga, gf, self.results.data.reactions.nua)
            self.results = self.strategy.solver()
            return self.results.j

        except AttributeError as e:
            print(f"Attribute error in current_energies: {e}")
            raise

        except Exception as e:
            print(f"Unexpected error in current_energies: {e}")
            raise

    def display_error_evolution(self, xk, convergence=0):
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


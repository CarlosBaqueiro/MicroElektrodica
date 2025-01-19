"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Grapher class

"""

import os
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from . import Calculator
from .tools import Tool


# for debugging
# import sys
# sys.exit()
# from .Tools import showme


class Grapher:

    def __init__(self, results):
        self.results = copy.deepcopy(results)
        self.data = self.results.data

        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions

        #self.graph_results(self.operation, self.species, self.results)

    def me(self):
        return r"$\mu$Elektrodica"

    def graph_results(self, operation, species, results):
        self.plot_fval(operation, species, results)
        self.plot_species(operation, species, results, "inSolution", login=True)
        self.plot_species(operation, species, results, "coverages", login=True)
        self.plot_species(operation, species, results, "all", login=True)
        self.grah_each_theta(operation, species, results, login=True)
        self.current(operation, results)

    def plot_species(self, operation, species, results, label, login=False):
        if label == "coverages":
            c_species = results.theta
            legends = species.adsorbed
            plt.ylabel(r"$\theta_i$")
        elif label == "inSolution":
            c_species = np.concatenate(
                [results.c_reactants, results.c_products], axis=1
            )
            legends = species.reactants + species.products
            plt.ylabel(r"$c_i\ mol/L$")
        elif label == "all":
            c_species = np.concatenate(
                [results.c_reactants, results.c_products, results.theta], axis=1
            )
            legends = species.reactants + species.products + species.adsorbed
            plt.ylabel(r"$c_i\ and\ \theta_i$")
        return self.plot(operation, c_species, legends, login=login)

    def grah_each_theta(self, operation, species, results, login=True):
        for i in range(len(species.adsorbed)):
            plt.ylabel(r"$\theta_i$")
            c_species = results.theta[:, i]
            legends = species.adsorbed[i]
            self.plot(operation, c_species, legends, login=login)

    def plot2(self, operation, c_species, legends, login=False):
        plt.plot(operation.potential, c_species, label=legends)
        plt.xlabel("Potential [V]")
        plt.grid(
            visible=True,
            which="both",
            axis="both",
            color="grey",
            linestyle="-",
            linewidth="0.2",
        )
        plt.minorticks_on()
        plt.legend(loc="lower right")
        if login:
            plt.yscale("log")
        plt.tight_layout()
        plt.show()

    def current(self, operation, results):
        plt.plot(results.j, operation.potential)
        plt.xlabel("Current density [A/cm2]")
        plt.ylabel("Potencial [V]")
        # plt.xscale('log')
        plt.tight_layout()
        plt.grid(
            visible=True,
            which="both",
            axis="both",
            color="grey",
            linestyle="-",
            linewidth="0.2",
        )
        plt.show()

    def plot_fval(self, operation, species, results):
        if operation.cstr == False:
            legends = Tool.format_latex_chemical(species.adsorbed)
            plt.ylabel(r"$\frac{\partial \theta_i}{\partial t}$")
        else:
            legends = Tool.format_latex_chemical(species.reactants + species.products + species.adsorbed)
            plt.ylabel(
                r"$\frac{\partial \c_i}{\partial t}\ and\ \frac{\partial \theta_i}{\partial t}$"
            )
        plt.plot(operation.potential, results.fval, linestyle="-.", label=legends)
        plt.xlabel("Potential [V]")
        ax = plt.gca()
        ax.yaxis.label.set_size(15)
        plt.yscale("log")
        plt.grid(
            True,
            which="both",
            axis="both",
            color="grey",
            linestyle="-",
            linewidth="0.2",
        )
        plt.minorticks_on()
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_experimental_points(self, ax, collection):
        for label, (marker, linestyle, data, color) in collection.items():
            ax.plot(
                data.potential,
                data.variable,
                label=label,
                marker=marker,
                linestyle=linestyle,
                color=color,
            )

    def plot_results(self, ax, collection, variable, subvariable=None):
        for calculator, (marker, linestyle, color) in collection.items():
            potential = calculator.operation.potential
            label = Tool.format_latex_chemical(calculator.data.species[subvariable])
            ax.plot(potential, calculator.results[variable],
                    label=label,
                    marker=marker,
                    linestyle=linestyle,
                    color=color,
                    )

    def plot_coverages(self, collection, figname):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.plot_results(ax, collection, "theta", "adsorbed")
        ax.set_title("Coverages")
        ax.set_xlabel("Potential [V]")
        ax.set_ylabel(r"$\theta_j$")
        ax.grid(True)
        fname = os.path.join(self.data.directory, figname)
        fig.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.close(fig)

    def plot_concentrations(self, collection, figname, role='Reactants'):
        global subvariable
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if role == 'Reactants':
            subvariable = 'reactants'
        elif role == 'Products':
            subvariable = 'products'
        self.plot_results(ax, collection, 'c_' + subvariable, subvariable)
        ax.set_title("Coverages")
        ax.set_xlabel("Potential [V]")
        ax.set_ylabel(r"$\theta_j$")
        ax.grid(True)
        fname = os.path.join(self.data.directory, figname)
        fig.savefig(fname, dpi=300, bbox_inches="tight", format="png")

"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Grapher class

"""

import matplotlib.pyplot as plt
import numpy as np

# for debugging
#import sys
# sys.exit()
#from .Tools import showme

class Grapher:
    """
    Provides methods for graphing chemical species and reactions based on the given data.

    The class includes functions to plot various graphs such as species concentration,
    reaction progress, and current density over potential. Results are visualized using
    matplotlib to analyze chemical reactions and species behavior under different conditions.

    :ivar operation: Defines the operation parameters, including potential and other settings.
    :type operation: dict
    :ivar species: Contains information about the different species involved in the reaction.
    :type species: dict
    :ivar reactions: Describes the chemical reactions present in the data.
    :type reactions: dict
    """
    def __init__(self, data, results):
        self.operation = data.parameters
        self.species = data.species
        self.reactions = data.reactions
        self.results = results.results

        self.graph_results(self.operation, self.species, self.results)

    def graph_results(self, operation, species, results):
        self.plot_fval(operation, species, results)
        self.plot_species(operation, species, results, 'inSolution', login=True)
        self.plot_species(operation, species, results, 'coverages', login=True)
        self.plot_species(operation, species, results, 'all', login=True)
        self.grah_each_theta(operation, species, results, login=True)
        self.current(operation, results)



    def plot_species(self, operation, species, results, label, login=False):
         if label == 'coverages':
            c_species = results.theta
            legends = species.adsorbed
            plt.ylabel(r'$\theta_i$')
         elif label == 'inSolution':
            c_species = np.concatenate([results.c_reactants, results.c_products], axis=1)
            legends = species.reactants + species.products
            plt.ylabel(r'$c_i\ mol/L$')
         elif label == 'all':
            c_species = np.concatenate([results.c_reactants, results.c_products, results.theta], axis=1)
            legends = species.reactants + species.products + species.adsorbed
            plt.ylabel(r'$c_i\ and\ \theta_i$')
         return self.plot(operation, c_species, legends, login=login)


    def grah_each_theta(self, operation, species, results, login=True):
        for i in range(len(species.adsorbed)):
            plt.ylabel(r'$\theta_i$')
            c_species = results.theta[:, i]
            legends = species.adsorbed[i]
            self.plot(operation, c_species, legends, login=login)

    def plot_fval(self, operation, species, results):
        if operation.cstr == False:
            legends = species.adsorbed
            plt.ylabel(r'$\frac{\partial \theta_i}{\partial t}$')
        else:
            legends = species.reactants + species.products + species.adsorbed
            plt.ylabel(r'$\frac{\partial \c_i}{\partial t}\ and\ \frac{\partial \theta_i}{\partial t}$')
        plt.plot(operation.potential, results.fval, linestyle='-.', label=legends)
        plt.xlabel('Potential [V]')
        ax = plt.gca()
        ax.yaxis.label.set_size(15)
        plt.yscale('log')
        plt.grid(True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot(self, operation, c_species, legends, login=False):
         plt.plot(operation.potential, c_species,  label=legends)
         plt.xlabel('Potential [V]')
         plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
         plt.minorticks_on()
         plt.legend(loc='lower right')
         if login: plt.yscale('log')
         plt.tight_layout()
         plt.show()

    def current(self, operation, results):
        plt.plot(results.j, operation.potential)
        plt.xlabel('Current density [A/cm2]')
        plt.ylabel('Potencial [V]')
        #plt.xscale('log')
        plt.tight_layout()
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.show()
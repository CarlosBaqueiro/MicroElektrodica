"""
M-elektrodica:
                Graphics functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""


import matplotlib.pyplot as plt
import numpy as np


class Grapher:
    def __init__(self, data, results):
        self.operation = data.parameters
        self.species = data.species
        self.reactions = data.reactions

        self.graph_results(self.operation, self.species, results)

    def graph_results(self, operation, species, results):
        self.plot_fval(operation, species, results)
        #self.plot_species(operation, species, results, 'inSolution', login=False)
        self.plot_species(operation, species, results, 'coverages', login=False, normal=False)
        #self.plot_species(operation, species, results, 'all', login=False)
        #self.grah_each_theta(operation, species, results, login=False)
        #self.current(operation, results)



    def plot_species(self, operation, species, results, label, login=False, normal = True):
         if label == 'coverages':
            c_species = results.theta
            legends = species.adsorbed
            plt.ylabel(r'$\theta_i$')
            if normal:
                sum_species = np.sum(c_species, axis=1)
                c_species = c_species/np.max(c_species, axis=1)[:, None]
                print(normal)
         elif label == 'inSolution':
            c_species = np.concatenate([results.c_reactants, results.c_products], axis=1)
            legends = species.reactants + species.products
            plt.ylabel(r'$c_i\ mol/L$')
         elif label == 'all':
            c_species = np.concatenate([results.c_reactants, results.c_products, results.theta], axis=1)
            legends = species.reactants + species.products + species.adsorbed
            plt.ylabel(r'$\frac{\partial \c_i}{\partial t}\ and\ \frac{\partial \theta_i}{\partial t}$')
         self.plot(operation, c_species, legends, login=False)
         self.plot(operation, c_species, legends, login=True)

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
         sum_species = np.sum(c_species, axis=1)
         plt.plot(operation.potential, sum_species,  label='sum')
         plt.xlabel('Potential [V]')
         plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
         plt.minorticks_on()
         plt.legend(loc='best')
         if login: plt.yscale('log')
         plt.tight_layout()
         plt.show()

    def current(self, operation, results):
        plt.plot(results.j, operation.potential)
        plt.xlabel('Current density [A/cm2]')
        plt.ylabel('Ptencial [V]')
        #plt.xscale('log')
        plt.tight_layout()
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.show()


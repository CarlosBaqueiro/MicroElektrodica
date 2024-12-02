"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Examples/Oxygen.py

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from Source import Collector, Calculator

# for debugging
#from .Tools import showme
#import sys
# sys.exit()

class Oxygen:
    def __init__(self):
        directory = os.path.join('Moore2013Oxygen')
        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions
        self.melek = Calculator(self.data)

        # Moore's Result
        o = np.array([
            [8.99231e-2, 1.64815e-2],
            [1.00532e-1, 1.76392e-2],
            [1.21005e-1, 1.82384e-2],
            [1.41249e-1, 2.45038e-2],
            [1.51291e-1, 3.03143e-2],
            [1.61011e-1, 3.65082e-2],
            [1.75140e-1, 4.52557e-2],
            [1.86489e-1, 5.38836e-2],
            [2.04009e-1, 7.09880e-2],
            [2.19247e-1, 9.40500e-2],
            [2.29825e-1, 1.07589e-1],
            [2.41768e-1, 1.33539e-1],
            [2.56270e-1, 1.65131e-1],
            [2.68212e-1, 1.92209e-1],
            [2.78066e-1, 2.12101e-1],
            [2.86127e-1, 2.32827e-1],
            [2.94384e-1, 2.50284e-1],
            [3.01823e-1, 2.67352e-1],
            [3.12572e-1, 2.92814e-1],
            [3.18543e-1, 3.11806e-1],
            [3.27647e-1, 3.33292e-1],
            [3.45847e-1, 3.37442e-1],
            [3.65931e-1, 3.35782e-1],
            [3.79738e-1, 3.34122e-1],
            [3.88524e-1, 3.12540e-1],
            [3.97311e-1, 2.95938e-1],
            [4.01076e-1, 2.84317e-1],
            [4.07352e-1, 2.71036e-1],
            [4.11118e-1, 2.57755e-1],
            [4.17394e-1, 2.41153e-1],
            [4.27435e-1, 2.17911e-1],
            [4.34967e-1, 2.02970e-1],
            [4.42498e-1, 1.83048e-1],
            [4.51284e-1, 1.68107e-1],
            [4.58816e-1, 1.44865e-1],
            [4.67602e-1, 1.28263e-1],
            [4.78899e-1, 1.08341e-1],
            [4.90196e-1, 9.00798e-2],
            [5.01492e-1, 7.59685e-2],
            [5.11534e-1, 6.35173e-2],
            [5.26596e-1, 4.85760e-2]
        ])
        oh = np.array([
            [1.32576e-1, 9.75229e-1],
            [1.11124e-1, 9.79091e-1],
            [8.72756e-2, 9.82411e-1],
            [1.55779e-1, 9.64623e-1],
            [1.77788e-1, 9.50407e-1],
            [1.98262e-1, 9.27841e-1],
            [2.07645e-1, 9.10917e-1],
            [2.26413e-1, 8.84591e-1],
            [2.38527e-1, 8.55409e-1],
            [2.46886e-1, 8.27049e-1],
            [2.56270e-1, 8.00347e-1],
            [2.68025e-1, 7.69912e-1],
            [2.76811e-1, 7.40859e-1],
            [2.84421e-1, 7.10273e-1],
            [2.88686e-1, 6.81126e-1],
            [2.94657e-1, 6.48030e-1],
            [3.01482e-1, 6.18069e-1],
            [3.07453e-1, 5.85694e-1],
            [3.12572e-1, 5.61154e-1],
            [3.18543e-1, 5.25801e-1],
            [3.22808e-1, 4.97782e-1],
            [3.29633e-1, 4.69200e-1],
            [3.35806e-1, 4.35391e-1],
            [3.39869e-1, 4.05829e-1],
            [3.45841e-1, 3.74989e-1],
            [3.53567e-1, 3.43019e-1],
            [3.57772e-1, 3.11710e-1],
            [3.62049e-1, 2.82847e-1],
            [3.68873e-1, 2.54264e-1],
            [3.73992e-1, 2.23237e-1],
            [3.79963e-1, 1.89201e-1],
            [3.87641e-1, 1.64003e-1],
            [3.99584e-1, 1.32129e-1],
            [4.08114e-1, 1.07589e-1],
            [4.20057e-1, 7.48694e-2],
            [4.30470e-1, 4.81214e-2],
            [4.50657e-1, 3.11444e-2],
            [4.71995e-1, 1.62030e-2],
            [4.93961e-1, 7.07221e-3],
            [5.17810e-1, 4.58198e-3],
            [5.41659e-1, 1.26168e-3],
            [5.64880e-1, 4.31606e-4],
            [5.88729e-1, 4.31606e-4],
            [6.12578e-1, 4.31606e-4],
            [6.36427e-1, 4.31606e-4]
        ])

        species = {
            'OH': ('s', '-', oh, 1),
            'O': ('o', '--', o, 0)
        }
        colors = ['#ffd9bf', '#ffb380', '#ff8c40', '#ff7f0e', '#ff6610']

        fname = os.path.join(directory, 'MooreOxygen_theta.png')
        for label, (mark, linest, data, idx) in species.items():
            plt.plot(data[:, 0] + 0.02, data[:, 1], marker=mark, linestyle='', color='#1f77b4')

        co = np.array([0.01, 0.1, 0.5, 1, 5])
        j = np.zeros((len(self.operation.potential), len(co)))
        pressure_legend = []
        for i in range(len(co)):
            self.data.species.c0_reactants = np.array([co[i], 1])
            self.melek = Calculator(self.data)
            plt.plot(self.operation.potential, self.melek.results.theta[:, 0], linestyle='--', color=colors[i])
            plt.plot(self.operation.potential, self.melek.results.theta[:, 1], linestyle='-', color=colors[i])
            j[:, i] = self.melek.results.j
            pressure_legend.append(Line2D([0], [0], linestyle='-', color=colors[i], label=fr'{co[i]} atm'))

        plt.xlabel('Overpotential [V]')
        plt.ylabel(r'$\theta_i$')
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        species_legend = [
            Line2D([0], [0], marker='o', markerfacecolor='#1f77b4', markeredgecolor='#1f77b4',
                   linestyle='--', color='#ff6610', label=fr'O*'),
            Line2D([0], [0], marker='s', markerfacecolor='#1f77b4', markeredgecolor='#1f77b4',
                   linestyle='-', color='#ff7f0e', label=fr'OH*')
        ]
        solutions_legend = [
            Line2D([0], [0], color='#1f77b4', lw=2, linestyle='-', label='Moore et. al. 2017'),
            Line2D([0], [0], color='#ff7f0e', lw=2, linestyle='-', label=r'$\mu$Elektrodica')
        ]

        first_legend = plt.legend(handles=species_legend, loc='lower right')
        second_legend = plt.legend(handles=pressure_legend, loc='center right')
        plt.gca().add_artist(first_legend)
        plt.gca().add_artist(second_legend)
        plt.legend(handles=solutions_legend, loc='upper right')
        plt.ylim(0, 1)
        plt.xlim(0.1, 0.9)
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
        plt.show()

        fname = os.path.join(directory, 'MooreOxygen_j.png')
        for i in range(len(co)):
            plt.plot(abs(j[:, i]), self.operation.potential, color=colors[i], label=f'{co[i]} atm')

        plt.xlabel(r'Current density [Acm$^{-2}$]')
        plt.ylabel('Overpotential [V]')
        plt.xscale('log')
        # plt.yscale('log')
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.legend(loc='lower right')
        plt.ylim(0.1, 0.9)
        plt.xlim(1e-7, 1e1)
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
        plt.show()

if __name__ == '__main__':
    Oxygen()
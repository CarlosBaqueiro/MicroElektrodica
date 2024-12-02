"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Examples/Hydrogen.py

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

class Hydrogen:
    def __init__(self):
        directory = os.path.join('Wang2007Hydrogen')

        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions

        self.melek_hb = Calculator(self.data)

        # Atop
        self.data.species.G_formation_ads = np.array([75]) * 1e-3
        self.data.reactions.Ga = np.array([196, 294, 48]) * 1e-3
        self.melek_atop = Calculator(self.data)

        # Wang's Results
        h_hb = np.array([
            [7.92752e-03, 6.55677e-01],
            [2.49151e-02, 6.68118e-01],
            [4.07701e-02, 6.65630e-01],
            [5.77576e-02, 6.65630e-01],
            [7.19139e-02, 6.54432e-01],
            [8.89015e-02, 6.38258e-01],
            [1.04190e-01, 6.19596e-01],
            [1.20045e-01, 6.02177e-01],
            [1.37033e-01, 5.74806e-01],
            [1.52888e-01, 5.46190e-01],
            [1.68743e-01, 5.13841e-01],
            [1.84032e-01, 4.80249e-01],
            [2.00453e-01, 4.41680e-01],
            [2.16874e-01, 4.01866e-01],
            [2.32729e-01, 3.54588e-01],
            [2.48018e-01, 3.07309e-01],
            [2.63873e-01, 2.62519e-01],
            [2.80861e-01, 2.17729e-01],
            [2.96149e-01, 1.72939e-01],
            [3.12005e-01, 1.36858e-01],
            [3.28992e-01, 1.03266e-01],
            [3.44281e-01, 7.83826e-02],
            [3.60136e-01, 5.84759e-02],
            [3.75991e-01, 4.23017e-02],
            [3.92412e-01, 3.60809e-02],
            [4.08267e-01, 2.36392e-02],
            [4.25821e-01, 1.74184e-02],
            [4.41676e-01, 1.21742e-02],
            [4.56399e-01, 8.70918e-03],
            [4.73952e-01, 4.97667e-03],
            [4.89807e-01, 3.24417e-03]
        ])
        h_atop = np.array([
            [3.70714e-4, 5.00000e-2],
            [9.63855e-3, 3.94882e-2],
            [1.96478e-2, 3.79668e-2],
            [2.96571e-2, 2.80775e-2],
            [3.96664e-2, 2.17151e-2],
            [4.96756e-2, 2.05394e-2],
            [6.96942e-2, 1.40387e-2],
            [8.97127e-2, 1.00277e-2],
            [1.29379e-1, 6.70816e-3],
            [1.69416e-1, 3.87275e-3],
            [2.49861e-1, 1.17566e-3]
        ])

        j_hb = np.array([
            [9.94083e+00, 1.12183e-03],
            [1.23077e+01, 1.71343e-03],
            [1.84615e+01, 2.51158e-03],
            [2.60355e+01, 3.73200e-03],
            [3.50296e+01, 5.17895e-03],
            [4.59172e+01, 7.28531e-03],
            [5.72781e+01, 9.44146e-03],
            [6.95858e+01, 1.24037e-02],
            [8.14201e+01, 1.67471e-02],
            [9.46746e+01, 2.20007e-02],
            [1.07929e+02, 2.66272e-02],
            [1.20237e+02, 3.40383e-02],
            [1.33491e+02, 4.35108e-02],
            [1.47219e+02, 5.41188e-02],
            [1.60947e+02, 6.73132e-02],
            [1.74201e+02, 8.37258e-02],
            [1.88402e+02, 1.01329e-01],
            [2.01657e+02, 1.24325e-01],
            [2.15858e+02, 1.50464e-01],
            [2.30059e+02, 1.79628e-01],
            [2.44734e+02, 2.08660e-01],
            [2.60355e+02, 2.39086e-01],
            [2.74556e+02, 2.66576e-01],
            [2.89704e+02, 3.01307e-01],
            [3.05325e+02, 3.26878e-01],
            [3.20947e+02, 3.54619e-01],
            [3.36095e+02, 3.79499e-01],
            [3.52189e+02, 4.00600e-01],
            [3.67811e+02, 4.22881e-01],
            [3.83432e+02, 4.34366e-01],
            [3.99527e+02, 4.46156e-01],
            [4.14675e+02, 4.58280e-1]
        ])
        j_atop = np.array([
            [1.06724e+0, 2.97346e-2],
            [2.13447e+0, 4.57539e-2],
            [2.98826e+0, 7.21587e-2],
            [5.33618e+0, 1.24812e-1],
            [8.75133e+0, 1.87382e-1],
            [1.45144e+1, 2.79592e-1],
            [2.00640e+1, 3.44694e-1],
            [2.92423e+1, 4.14617e-1],
            [4.20491e+1, 4.71842e-1],
            [5.78442e+1, 4.92622e-1],
            [7.23586e+1, 5.17492e-1],
            [8.62327e+1, 5.27139e-1],
            [9.90395e+1, 5.36966e-1],
            [1.10566e+2, 5.60613e-1],
            [1.22732e+2, 5.78139e-1],
            [1.35752e+2, 6.14851e-1],
            [1.47919e+2, 6.45893e-1],
            [1.64781e+2, 7.21587e-1],
            [1.89755e+2, 8.78720e-1],
            [2.06190e+2, 1.03126e+0],
            [2.27748e+2, 1.32738e+0],
            [2.48026e+2, 1.75112e+0],
            [2.61046e+2, 2.13244e+0],
            [2.72572e+2, 2.53363e+0],
            [2.83671e+2, 2.97346e+0],
            [2.98613e+2, 4.04531e+0],
            [3.16969e+2, 5.33670e+0],
            [3.33618e+2, 7.30527e+0]
        ])

        species = {
            'H/B': ('o', ':', h_hb),
            'Atop': ('x', '-', h_atop)
        }

        fname = os.path.join(directory, 'WangHydrogen_theta.png')
        plt.plot(self.operation.potential, self.melek_hb.results.theta, linestyle=':', color='#ff7f0e')
        plt.plot(self.operation.potential, self.melek_atop.results.theta, linestyle='-', color='#ff7f0e')
        for label, (mark, linest, data) in species.items():
            plt.plot(data[:, 0], data[:, 1], marker=mark, linestyle='', color='#1f77b4')
        plt.xlabel('Overpotential [V]')
        plt.ylabel(r'$\theta_i$')
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()

        species_legend = [
            Line2D([0], [0], marker='o', markerfacecolor='#1f77b4', markeredgecolor='#1f77b4',
                   linestyle=':', color='#ff7f0e', label=fr'H/B'),
            Line2D([0], [0], marker='x', markerfacecolor='#1f77b4', markeredgecolor='#1f77b4',
                   linestyle='-', color='#ff7f0e', label=fr'Atop')
        ]
        solutions_legend = [
            Line2D([0], [0], color='#1f77b4', lw=2, linestyle='-', label='Wang et. al. 2007'),
            Line2D([0], [0], color='#ff7f0e', lw=2, linestyle='-', label=r'$\mu$Elektrodica')
        ]
        first_legend = plt.legend(handles=species_legend, loc='center left')
        plt.gca().add_artist(first_legend)
        plt.legend(handles=solutions_legend, loc='center right')
        plt.tight_layout()
        plt.xlim(0, 0.5)
        plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
        plt.show()

        currents = {
            'H/B': ('o', ':', j_hb),
            'Atop': ('x', '-', j_atop)
        }

        fname = os.path.join(directory, 'WangHydrogen_j.png')
        plt.plot(self.operation.potential, self.melek_hb.results.j, linestyle=':', color='#ff7f0e')
        plt.plot(self.operation.potential, self.melek_atop.results.j, linestyle='-', color='#ff7f0e')
        for label, (mark, linest, data) in currents.items():
            plt.plot(data[:, 0] * 1e-3, data[:, 1], marker=mark, linestyle='', color='#1f77b4')
        plt.ylabel(r'Current density [Acm$^{-2}$]')
        plt.xlabel('Overpotential [V]')
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()
        plt.yscale('log')
        plt.ylim(1e-3, 10)
        plt.xlim(0, 0.4)
        first_legend = plt.legend(handles=species_legend, loc='upper left')
        plt.gca().add_artist(first_legend)
        plt.legend(handles=solutions_legend, loc='lower right')
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
        plt.show()

if __name__ == '__main__':
    Hydrogen()
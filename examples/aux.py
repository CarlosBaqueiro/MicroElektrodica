"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        examples/hydrogen.py

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.lines import Line2D

from melektrodica import Collector, Fitter, Calculator, CoordinatorReactions, Kpynetic, Grapher
from melektrodica.constants import F, k_B, h


# for debugging
# from .Tools import showme
# import sys
# sys.exit()


class Hydrogen:
    def __init__(self):
        directory = os.path.join("Wang2007Hydrogen")

        self.data = Collector(directory)
        _data = copy.deepcopy(self.data)
        _melek = Calculator(_data)

        self.data.species.g_formation_ads = self.data.species.g_formation_ads
        self.data.reactions.ga = self.data.reactions.ga

        print('\ng0:')
        print('ga: ', self.data.reactions.ga)
        print('g_formation: ', self.data.species.g_formation_ads)

        j_hb = np.array(
            [
                [9.94083e00, 1.12183e-03],
                [1.23077e01, 1.71343e-03],
                [1.84615e01, 2.51158e-03],
                [2.60355e01, 3.73200e-03],
                [3.50296e01, 5.17895e-03],
                [4.59172e01, 7.28531e-03],
                [5.72781e01, 9.44146e-03],
                [6.95858e01, 1.24037e-02],
                [8.14201e01, 1.67471e-02],
                [9.46746e01, 2.20007e-02],
                [1.07929e02, 2.66272e-02],
                [1.20237e02, 3.40383e-02],
                [1.33491e02, 4.35108e-02],
                [1.47219e02, 5.41188e-02],
                [1.60947e02, 6.73132e-02],
                [1.74201e02, 8.37258e-02],
                [1.88402e02, 1.01329e-01],
                [2.01657e02, 1.24325e-01],
                [2.15858e02, 1.50464e-01],
                [2.30059e02, 1.79628e-01],
                [2.44734e02, 2.08660e-01],
                [2.60355e02, 2.39086e-01],
                [2.74556e02, 2.66576e-01],
                [2.89704e02, 3.01307e-01],
                [3.05325e02, 3.26878e-01],
                [3.20947e02, 3.54619e-01],
                [3.36095e02, 3.79499e-01],
                [3.52189e02, 4.00600e-01],
                [3.67811e02, 4.22881e-01],
                [3.83432e02, 4.34366e-01],
                [3.99527e02, 4.46156e-01],
                [4.14675e02, 4.58280e-1],
            ]
        )

        vdata = j_hb[:, 0] * 1e-3
        jdata = j_hb[:, 1]

        F = Fitter(self.data, vdata, jdata)
        self.data.reactions.ga = F.ga_fit
        self.data.species.g_formation_ads = F.gf_fit
        print('\ng_fit:')
        print('ga: ', self.data.reactions.ga)
        print('g_formation: ', self.data.species.g_formation_ads)

        self.nmelek = Calculator(self.data)

        plt.plot(vdata, jdata, 'o', label='Experimental Data')
        #plt.plot(_data.parameters.potential, abs(_melek.results.j), label='Original')
        plt.plot(self.data.parameters.potential, abs(self.nmelek.results.j), label='Fit')
        plt.xlabel('Potential (V)')
        plt.ylabel('J (A/cm$^2$)')
        plt.yscale('log')
        plt.ylim(1e-3, 10)
        plt.xlim(0, 0.5)
        plt.legend()
        plt.show()


class Hydrogen2:
    def __init__(self):
        directory = os.path.join("Wang2007Hydrogen")
        self.data = Collector(directory)
        self.data.species.g_formation_ads = np.array([75]) * 1e-3
        self.data.reactions.ga = np.array([196, 294, 48]) * 1e-3

        source = self.data.species.reactants[0]
        target = self.data.species.products[0]
        eta = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        fname = "Wang2007Hydrogen"
        self.Kpy = Kpynetic(self.data)
        coo = CoordinatorReactions(self.Kpy)
        Coord = coo.plot_rxn_coords_potential(source, target, eta, fname)
        print(k_B * self.data.parameters.temperature / h)
        print(k_B * self.data.parameters.temperature * F * 2 / h)


class Oxygen:
    def __init__(self):
        directory = os.path.join("Moore2013Oxygen")
        self.data = Collector(directory)
        self.data.species.c0_reactants = np.array([0.21, 1])
        self.data.parameters.temperature = 50 + 273.15

        _data = copy.deepcopy(self.data)
        _Kpy = Kpynetic(_data)
        _melek = Calculator(_Kpy)

        print('g0:')
        print('ga: ', self.data.reactions.ga)
        print('g_formation: ', self.data.species.g_formation_ads)

        p021 = np.array(
            [
                [1.12533e-5, 9.62268e-1],
                [1.40802e-5, 9.55855e-1],
                [1.69398e-5, 9.50836e-1],
                [2.00405e-5, 9.45818e-1],
                [2.35762e-5, 9.41078e-1],
                [2.77357e-5, 9.36059e-1],
                [3.28124e-5, 9.31320e-1],
                [3.83857e-5, 9.26580e-1],
                [4.59235e-5, 9.21283e-1],
                [5.52502e-5, 9.15985e-1],
                [6.72201e-5, 9.10409e-1],
                [8.22428e-5, 9.04275e-1],
                [1.08834e-4, 8.95353e-1],
                [1.42417e-4, 8.86710e-1],
                [1.91657e-4, 8.77230e-1],
                [2.69747e-4, 8.66078e-1],
                [3.49050e-4, 8.57993e-1],
                [4.31870e-4, 8.50465e-1],
                [6.04437e-4, 8.36245e-1],
                [7.23130e-4, 8.28439e-1],
                [9.14987e-4, 8.18401e-1],
                [1.15775e-3, 8.08364e-1],
                [1.35440e-3, 8.01952e-1],
                [1.63862e-3, 7.90242e-1],
                [2.02742e-3, 7.77138e-1],
                [2.34535e-3, 7.67937e-1],
                [2.83753e-3, 7.56506e-1],
                [3.47168e-3, 7.43123e-1],
                [4.10713e-3, 7.31691e-1],
                [5.08163e-3, 7.16636e-1],
                [6.28736e-3, 7.02416e-1],
                [7.52201e-3, 6.90428e-1],
                [8.99911e-3, 6.79275e-1],
                [1.07663e-2, 6.67565e-1],
                [1.23158e-2, 6.59201e-1],
                [1.33208e-2, 6.54182e-1],
                [1.49003e-2, 6.46933e-1],
                [1.70449e-2, 6.39405e-1],
                [2.18102e-2, 6.24349e-1],
                [2.69852e-2, 6.10688e-1],
            ]
        )
        vdata = 1.18 - p021[:, 1]
        jdata = p021[:, 0]

        F = Fitter(self.data, vdata, jdata)
        self.data.reactions.ga = F.ga_fit
        self.data.species.g_formation_ads = F.gf_fit
        print('\ng_fit:')
        print('ga: ', self.data.reactions.ga)
        print('g_formation: ', self.data.species.g_formation_ads)
        self.data.parameters.potential = _data.parameters.potential
        self.nmelek = Calculator(self.data)
        G = Grapher(self.nmelek.results)
        list = {self.melek, ('', '-', 'tab:orange')}
        G.plot_results(list, 'prueba')
        plt.plot(vdata, jdata, 'o', label='Experimental Data')
        plt.plot(_data.parameters.potential, abs(_melek.results.j), label='Original')
        plt.plot(self.data.parameters.potential, abs(self.nmelek.results.j), label='Fit')
        plt.xlabel('Potential (V)')
        plt.ylabel('J (A/cm$^2$)')
        plt.yscale('log')
        plt.legend()
        plt.show()

        plt.plot(_data.parameters.potential, _melek.results.theta, label='Original')
        plt.plot(self.data.parameters.potential, self.nmelek.results.theta, label='Fit')
        plt.xlabel('Potential (V)')
        plt.ylabel('fc')
        plt.yscale('log')
        plt.legend()
        plt.show()

        plt.plot(_data.parameters.potential, _melek.results.theta, label='Original')
        plt.plot(self.data.parameters.potential, self.nmelek.results.theta, label='Fit')
        plt.xlabel('Potential (V)')
        plt.ylabel('fc')
        plt.legend()
        plt.show()


class Oxygen2:
    def __init__(self):
        directory = os.path.join("Moore2013Oxygen")
        self.data = Collector(directory)
        #self.data.reactions.ga = np.array([258, 459, 502, 455]) * 1e-3
        #self.data.species.g_formation_ads = np.array([-477, -120]) * 1e-3
        coo = CoordinatorReactions(self.data)
        source = self.data.species.reactants[0]
        target = self.data.species.products[0]
        eta = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        fname = "Moore2013Oxygen"
        Coord = coo.plot_rxn_coords_potential(source, target, eta, fname)


if __name__ == "__main__":
    #Hydrogen()
    #Hydrogen2()
    Oxygen()
    #Oxygen2()

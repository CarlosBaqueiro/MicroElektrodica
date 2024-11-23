import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
#sys.exit()


from .Kpynetic import *
from .Collector import *
from .Calculator import *
from .Constants import k_B
from .Tools import showme


class Validation:
    def __init__(self):
        self.hydrogen = Hydrogen()
        self.oxygen = Oxygen()
        self.ethanol = Ethanol()
    
class Hydrogen:
    def __init__(self):
        directory = os.path.join('Examples', 'Wang2007Hydrogen')

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


class Oxygen:
    def __init__(self):
        directory = os.path.join('Examples', 'Moore2013Oxygen')
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
        
        
class Ethanol:
    def __init__(self):
        directory = os.path.join('Examples', 'SanchezMonreal2017Ethanol')
        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions
        self.melek = Calculator(self.data)
        
        # CSTR mode
        self.data.parameters.cstr = True
        self.data.parameters.Fv = 3e-4
        self.data.parameters.Ac = 1

        self.data.reactions.nux = self.data.reactions.nuc
        self.melek_cstr = Calculator(self.data)

        ch3choh =np.array([
            [1.04077e-1, 1.99737e-4],
            [1.11513e-1, 2.08246e-4],
            [1.18949e-1, 2.08246e-4],
            [1.26385e-1, 2.08246e-4],
            [1.33821e-1, 1.98353e-4],
            [1.41257e-1, 1.87620e-4],
            [1.48693e-1, 1.73202e-4],
            [1.56129e-1, 1.58232e-4],
            [1.63565e-1, 1.41572e-4],
            [1.71001e-1, 1.26666e-4],
            [1.78437e-1, 1.11377e-4],
            [1.85873e-1, 9.82741e-5],
            [1.93309e-1, 8.58130e-5],
            [2.00745e-1, 7.49320e-5],
            [2.08181e-1, 6.52036e-5],
            [2.15617e-1, 5.67382e-5],
            [2.23053e-1, 4.86900e-5],
            [2.30489e-1, 4.23686e-5],
            [2.37925e-1, 3.63587e-5],
            [2.45361e-1, 3.13100e-5],
            [2.52797e-1, 2.68688e-5],
            [2.60233e-1, 2.32184e-5],
            [2.67669e-1, 1.99249e-5],
            [2.75105e-1, 1.70986e-5],
            [2.82541e-1, 1.46223e-5],
            [2.89977e-1, 1.25046e-5],
            [2.97413e-1, 1.06565e-5],
            [3.04849e-1, 9.17676e-6],
            [3.12285e-1, 7.87506e-6],
            [3.19721e-1, 6.75800e-6],
            [3.27157e-1, 5.79940e-6],
            [3.34593e-1, 4.94228e-6],
            [3.42029e-1, 4.24123e-6],
            [3.49465e-1, 3.58936e-6],
            [3.56901e-1, 3.03768e-6],
            [3.64337e-1, 2.56187e-6],
            [3.71773e-1, 2.13075e-6],
            [3.79209e-1, 1.75990e-6],
            [3.86645e-1, 1.44352e-6],
            [3.94081e-1, 1.15958e-6],
            [4.01517e-1, 9.31488e-7],
            [4.08953e-1, 7.32815e-7],
            [4.16389e-1, 5.60702e-7],
            [4.23825e-1, 4.26040e-7],
            [4.31261e-1, 3.20360e-7],
            [4.38697e-1, 2.32664e-7],
            [4.46133e-1, 1.67220e-7]
        ])
        ch3co=np.array([
            [1.02680e-1, 1.27841e-2],
            [1.11333e-1, 1.59723e-2],
            [1.20885e-1, 1.98203e-2],
            [1.32010e-1, 2.46791e-2],
            [1.45158e-1, 3.09649e-2],
            [1.55272e-1, 3.61751e-2],
            [1.67184e-1, 4.28041e-2],
            [1.82579e-1, 5.20882e-2],
            [1.94828e-1, 6.05430e-2],
            [2.08538e-1, 7.12730e-2],
            [2.20843e-1, 8.17926e-2],
            [2.36115e-1, 9.68796e-2],
            [2.49331e-1, 1.11918e-1],
            [2.65243e-1, 1.32224e-1],
            [2.78324e-1, 1.51895e-1],
            [2.93157e-1, 1.76280e-1],
            [3.08261e-1, 2.05416e-1],
            [3.22150e-1, 2.34775e-1],
            [3.37321e-1, 2.71081e-1],
            [3.50806e-1, 3.06369e-1],
            [3.65303e-1, 3.47134e-1],
            [3.82159e-1, 3.99386e-1],
            [3.97262e-1, 4.47709e-1],
            [4.12163e-1, 4.97293e-1],
            [4.28683e-1, 5.50682e-1],
            [4.46618e-1, 6.05464e-1]
        ])
        co = np.array([
            [1.05040e-1, 9.84787e-1],
            [1.16165e-1, 9.81029e-1],
            [1.27628e-1, 9.74796e-1],
            [1.39090e-1, 9.69220e-1],
            [1.50552e-1, 9.63676e-1],
            [1.62015e-1, 9.55116e-1],
            [1.73477e-1, 9.50864e-1],
            [1.85276e-1, 9.44355e-1],
            [1.97076e-1, 9.35833e-1],
            [2.08538e-1, 9.26338e-1],
            [2.20001e-1, 9.16940e-1],
            [2.31463e-1, 9.06480e-1],
            [2.42925e-1, 8.93859e-1],
            [2.54388e-1, 8.79729e-1],
            [2.65850e-1, 8.65271e-1],
            [2.77312e-1, 8.48883e-1],
            [2.88775e-1, 8.30686e-1],
            [3.00237e-1, 8.09260e-1],
            [3.11699e-1, 7.87382e-1],
            [3.23162e-1, 7.63658e-1],
            [3.33950e-1, 7.37352e-1],
            [3.45412e-1, 7.08783e-1],
            [3.56874e-1, 6.77856e-1],
            [3.67663e-1, 6.44981e-1],
            [3.78451e-1, 6.12138e-1],
            [3.89239e-1, 5.76909e-1],
            [4.00027e-1, 5.39564e-1],
            [4.10815e-1, 5.02390e-1],
            [4.20929e-1, 4.65102e-1],
            [4.31043e-1, 4.38392e-1],
            [4.41156e-1, 4.01822e-1],
            [4.49585e-1, 3.71279e-1]
        ])
        ch3=np.array([
            [1.04842e-1, 1.35525e-4],
            [1.12263e-1, 7.96087e-5],
            [1.19685e-1, 4.62775e-5],
            [1.27107e-1, 2.69954e-5],
            [1.34528e-1, 1.59680e-5],
            [1.41950e-1, 9.54424e-6],
            [1.49372e-1, 5.78460e-6],
            [1.56793e-1, 3.50594e-6],
            [1.64215e-1, 2.12489e-6],
            [1.71637e-1, 1.30590e-6],
            [1.79058e-1, 8.05364e-7],
            [1.86480e-1, 4.93236e-7],
            [1.93902e-1, 3.08445e-7],
            [2.01323e-1, 1.91550e-7],
            [2.08745e-1, 1.18956e-7],
            [2.16167e-1, 7.46480e-8],
            [2.23588e-1, 4.63577e-8],
            [2.31010e-1, 2.89898e-8],
            [2.38432e-1, 1.81288e-8],
            [2.45853e-1, 1.12192e-8],
            [2.53275e-1, 6.96731e-9],
            [2.60696e-1, 4.29683e-9],
            [2.68118e-1, 2.65914e-9],
            [2.75540e-1, 1.64564e-9],
            [2.82961e-1, 1.01488e-9],
            [2.90383e-1, 6.28072e-10],
            [2.97805e-1, 3.85995e-10],
            [3.05226e-1, 2.38048e-10],
            [3.12648e-1, 1.46298e-10],
            [3.20070e-1, 8.95983e-11],
            [3.27491e-1, 5.46830e-11],
            [3.34913e-1, 3.34900e-11],
            [3.42335e-1, 2.02272e-11],
            [3.49756e-1, 1.23021e-11],
            [3.57178e-1, 7.45607e-12],
            [3.64600e-1, 4.47210e-12],
            [3.72021e-1, 2.69168e-12],
            [3.79443e-1, 1.60884e-12],
            [3.86865e-1, 9.61623e-13],
            [3.94286e-1, 5.70789e-13],
            [4.01708e-1, 3.37626e-13],
            [4.09130e-1, 1.99708e-13],
            [4.16551e-1, 1.16903e-13],
            [4.23973e-1, 6.81940e-14],
            [4.31394e-1, 3.96420e-14],
            [4.38816e-1, 2.28847e-14],
            [4.46238e-1, 1.32570e-14],
            [4.50623e-1, 1.00027e-14]
        ])
        oh =np.array([
            [1.04703e-1, 1.77543e-1],
            [1.11783e-1, 1.97609e-1],
            [1.21447e-1, 2.22574e-1],
            [1.27965e-1, 2.37424e-1],
            [1.38079e-1, 2.58043e-1],
            [1.46170e-1, 2.72466e-1],
            [1.54261e-1, 2.83810e-1],
            [1.61003e-1, 2.92627e-1],
            [1.72466e-1, 3.04034e-1],
            [1.85276e-1, 3.14280e-1],
            [1.98087e-1, 3.22395e-1],
            [2.09550e-1, 3.29038e-1],
            [2.19663e-1, 3.33542e-1],
            [2.32474e-1, 3.39260e-1],
            [2.41240e-1, 3.42737e-1],
            [2.52702e-1, 3.48020e-1],
            [2.64164e-1, 3.53385e-1],
            [2.72930e-1, 3.57007e-1],
            [2.85740e-1, 3.63126e-1],
            [2.95854e-1, 3.68097e-1],
            [3.08665e-1, 3.75682e-1],
            [3.18779e-1, 3.81473e-1],
            [3.30241e-1, 3.89333e-1],
            [3.41704e-1, 3.97355e-1],
            [3.51818e-1, 4.05543e-1],
            [3.65977e-1, 4.18141e-1],
            [3.78788e-1, 4.31131e-1],
            [3.90250e-1, 4.44524e-1],
            [4.04409e-1, 4.63032e-1],
            [4.17220e-1, 4.79857e-1],
            [4.22614e-1, 4.89744e-1],
            [4.34077e-1, 5.04959e-1],
            [4.46888e-1, 5.28671e-1]
        ])
        
        species = {
            'CH$_3$CHOH': ('tab:orange', ch3choh, 0),
            'CH$_3$CO': ('tab:blue', ch3co, 1),
            'CO': ('tab:green', co, 2),
            'CH$_3$': ('tab:red', ch3, 3),
            'OH': ('tab:purple', oh, 4)
        }

        fname = os.path.join(directory, 'SanchezEthanol_concentrations.png')
        plt.plot(self.operation.potential, self.melek_cstr.results.c_reactants[:,0], label=r'CH$_3$CH$_2$OH')
        plt.plot(self.operation.potential, self.melek_cstr.results.c_products[:,0], label=r'CH$_3$CHO')
        plt.xlabel('Overpotential [V]')
        plt.ylabel(r'$c_i$')
        plt.xlim(0, 1)
        #plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()
        plt.legend()
        plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
        plt.show()

        fname = os.path.join(directory, 'SanchezEthanol_theta.png')
        fname_zoom = os.path.join(directory, 'SanchezEthanol_theta_zoom.png')
        for lim in {True, False}:
            for label, (colors, data, idx) in species.items():
                plt.plot(data[:, 0], data[:, 1], linestyle='-', color=colors)
                plt.plot(self.operation.potential, self.melek.results.theta[:, idx], linestyle='--', color=colors)
                plt.plot(self.operation.potential, self.melek_cstr.results.theta[:, idx], linestyle=':', color=colors)

            plt.xlabel('Overpotential [V]')
            plt.ylabel(r'$\theta_i$')
            #plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
            plt.yscale('log')
            plt.minorticks_on()
            plt.tight_layout()

            species_legend = [
                Line2D([0], [0], color=colors, lw=2, label=fr'{specie}*') for specie, (colors, _, _) in species.items()
            ]
            solutions_legend = [
                Line2D([0], [0], color='grey', lw=2, linestyle='-', label='Sanchez-Monreal et. al. 2017'),
                Line2D([0], [0], color='grey', lw=2, linestyle='--', label=r'$\mu$Elektrodica'),
                Line2D([0], [0], color='grey', lw=2, linestyle=':', label=r'$\mu$Elektrodica CSTR')
            ]
            plt.xlim(0, 1)
            plt.ylim(1e-15, 1.5)
            if lim:
                plt.ylim(1e-3, 1.1)
                plt.legend(handles=solutions_legend, loc='lower center')
                plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
                plt.savefig(fname_zoom, dpi=300, bbox_inches='tight', format='png')
            else:
                first_legend = plt.legend(handles=species_legend, loc='lower left')
                plt.gca().add_artist(first_legend)
                #plt.legend(handles=solutions_legend, loc='lower right')
                plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
                plt.savefig(fname, dpi=300, bbox_inches='tight', format='png')
            plt.show()

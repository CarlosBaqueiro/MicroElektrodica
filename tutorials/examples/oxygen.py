"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        examples/oxygen.py

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from melektrodica import Collector, Calculator, Fitter


# for debugging
# from .Tools import showme
# import sys
# sys.exit()


class Oxygen:
    def __init__(self):
        directory = os.path.join("Moore2013Oxygen")
        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions
        self.melek = Calculator(self.data)

        # Moore's Result
        o = np.array(
            [
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
                [5.26596e-1, 4.85760e-2],
            ]
        )
        oh = np.array(
            [
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
                [6.36427e-1, 4.31606e-4],
            ]
        )

        species = {"OH": ("s", "-", oh, 1), "O": ("o", "--", o, 0)}
        colors = ["#ffd9bf", "#ffb380", "#ff8c40", "#ff7f0e", "#ff6610"]

        fname = os.path.join(directory, "MooreOxygen_theta.png")
        for label, (mark, linest, data, idx) in species.items():
            plt.plot(
                data[:, 0] + 0.02,
                data[:, 1],
                marker=mark,
                linestyle="",
                color="#1f77b4",
            )

        co = np.array([0.01, 0.1, 0.53, 1, 5])
        j = np.zeros((len(self.operation.potential), len(co)))
        pressure_legend = []
        for i in range(len(co)):
            self.data.species.c0_reactants = np.array([co[i], 1])
            self.melek = Calculator(self.data)
            plt.plot(
                self.operation.potential,
                self.melek.results.theta[:, 0],
                linestyle="--",
                color=colors[i],
            )
            plt.plot(
                self.operation.potential,
                self.melek.results.theta[:, 1],
                linestyle="-",
                color=colors[i],
            )
            j[:, i] = self.melek.results.j
            pressure_legend.append(
                Line2D([0], [0], linestyle="-", color=colors[i], label=rf"{co[i]} atm")
            )

        plt.xlabel("Overpotential [V]")
        plt.ylabel(r"$\theta_j$")
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        species_legend = [
            Line2D(
                [0],
                [0],
                marker="o",
                markerfacecolor="#1f77b4",
                markeredgecolor="#1f77b4",
                linestyle="--",
                color="#ff6610",
                label=rf"O*",
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                markerfacecolor="#1f77b4",
                markeredgecolor="#1f77b4",
                linestyle="-",
                color="#ff7f0e",
                label=rf"OH*",
            ),
        ]
        solutions_legend = [
            Line2D(
                [0],
                [0],
                color="#1f77b4",
                lw=2,
                linestyle="-",
                label="Moore et. al. 2017",
            ),
            Line2D(
                [0],
                [0],
                color="#ff7f0e",
                lw=2,
                linestyle="-",
                label=r"$\mu$Elektrodica",
            ),
        ]

        first_legend = plt.legend(handles=species_legend, loc="lower right")
        second_legend = plt.legend(handles=pressure_legend, loc="center right")
        plt.gca().add_artist(first_legend)
        plt.gca().add_artist(second_legend)
        plt.legend(handles=solutions_legend, loc="upper right")
        plt.ylim(0, 1)
        plt.xlim(0.1, 0.9)
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

        fname = os.path.join(directory, "MooreOxygen_j.png")
        for i in range(len(co)):
            plt.plot(
                self.operation.potential,
                abs(j[:, i]),
                color=colors[i],
                label=f"{co[i]} atm",
            )
        plt.ylabel(r"Current density [Acm$^{-2}$]")
        plt.xlabel("Overpotential [V]")
        plt.yscale("log")
        # plt.yscale('log')
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.legend(loc="lower right")
        plt.xlim(0.1, 0.9)
        plt.ylim(1e-7, 1e1)
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

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
        p105 = np.array(
            [
                [1.19018e-5, 9.93216e-1],
                [1.48085e-5, 9.87082e-1],
                [1.82197e-5, 9.80948e-1],
                [2.24168e-5, 9.75093e-1],
                [2.75807e-5, 9.68959e-1],
                [3.39342e-5, 9.62825e-1],
                [4.17512e-5, 9.56691e-1],
                [5.13690e-5, 9.50279e-1],
                [6.32023e-5, 9.44145e-1],
                [7.77616e-5, 9.37454e-1],
                [9.56747e-5, 9.31320e-1],
                [1.17057e-4, 9.24907e-1],
                [1.44022e-4, 9.18216e-1],
                [1.76208e-4, 9.11803e-1],
                [2.16799e-4, 9.05390e-1],
                [2.65251e-4, 8.99257e-1],
                [3.22717e-4, 8.92007e-1],
                [3.92634e-4, 8.84758e-1],
                [4.77697e-4, 8.77509e-1],
                [5.81190e-4, 8.70818e-1],
                [7.07104e-4, 8.63290e-1],
                [8.60298e-4, 8.56599e-1],
                [1.03502e-3, 8.49071e-1],
                [1.24522e-3, 8.40985e-1],
                [1.48142e-3, 8.32621e-1],
                [1.77233e-3, 8.24257e-1],
                [2.12036e-3, 8.16450e-1],
                [2.53673e-3, 8.07807e-1],
                [3.03487e-3, 8.00000e-1],
                [3.53052e-3, 7.90520e-1],
                [4.13020e-3, 7.81320e-1],
                [4.85889e-3, 7.72119e-1],
                [5.65244e-3, 7.62918e-1],
                [6.64969e-3, 7.53996e-1],
                [7.77918e-3, 7.44517e-1],
                [8.99911e-3, 7.35595e-1],
                [1.05868e-2, 7.25836e-1],
                [1.22470e-2, 7.16357e-1],
                [1.43273e-2, 7.07156e-1],
                [1.66672e-2, 6.97677e-1],
                [1.94982e-2, 6.88476e-1],
                [2.25559e-2, 6.79275e-1],
                [2.65354e-2, 6.70074e-1],
                [3.08691e-2, 6.60874e-1],
                [3.61124e-2, 6.51394e-1],
                [4.22463e-2, 6.42472e-1],
                [4.94221e-2, 6.33271e-1],
                [5.78167e-2, 6.24071e-1],
                [6.80172e-2, 6.15149e-1],
                [7.95703e-2, 6.05669e-1],
            ]
        )
        p188 = np.array(
            [
                [1.48917e-5, 9.98513e-1],
                [1.98171e-5, 9.89870e-1],
                [2.62244e-5, 9.81784e-1],
                [3.43166e-5, 9.73141e-1],
                [4.54118e-5, 9.65335e-1],
                [5.97586e-5, 9.57249e-1],
                [7.86379e-5, 9.48606e-1],
                [1.03482e-4, 9.39963e-1],
                [1.36174e-4, 9.31320e-1],
                [1.77198e-4, 9.22398e-1],
                [2.33180e-4, 9.13476e-1],
                [3.05133e-4, 9.04554e-1],
                [3.94840e-4, 8.95911e-1],
                [5.13790e-4, 8.85874e-1],
                [6.61125e-4, 8.75836e-1],
                [8.55491e-4, 8.66636e-1],
                [1.10081e-3, 8.56599e-1],
                [1.41648e-3, 8.46561e-1],
                [1.80237e-3, 8.35409e-1],
                [2.41199e-3, 8.21468e-1],
                [3.06907e-3, 8.10316e-1],
                [3.88334e-3, 7.99721e-1],
                [4.80474e-3, 7.87454e-1],
                [5.97818e-3, 7.75465e-1],
                [7.31421e-3, 7.63197e-1],
                [9.10051e-3, 7.51487e-1],
                [1.13231e-2, 7.39498e-1],
                [1.39315e-2, 7.27230e-1],
                [1.72370e-2, 7.14963e-1],
                [2.13269e-2, 7.02695e-1],
                [2.60931e-2, 6.90428e-1],
                [3.22842e-2, 6.78160e-1],
                [3.77679e-2, 6.69238e-1],
                [4.64681e-2, 6.56970e-1],
                [5.68530e-2, 6.45260e-1],
                [7.11353e-2, 6.32714e-1],
                [8.75220e-2, 6.19888e-1],
            ]
        )

        po2 = {
            "p021": ("s", "-", p021),
            #"p105": ("o", "--", p105),
            #"p188": ("x", "--", p188),
        }
        fname = os.path.join(directory, "MooreOxygen_j_50.png")
        for label, (mark, linest, data) in po2.items():
            plt.plot(
                1.18 - data[:, 1],
                data[:, 0],
                marker=mark,
                linestyle="",
                color="#1f77b4",
            )

        co = np.array([0.21, 1.05, 1.88])
        pressure_legend = []
        self.operation.temperature = 80 + 273.15
        for i in range(len(co)):
            self.data.species.c0_reactants = np.array([co[i], 1])
            self.melek = Calculator(self.data)
            plt.plot(
                self.operation.potential,
                abs(self.melek.results.j),
                linestyle="-",
                color=colors[i + 1],
            )
            pressure_legend.append(
                Line2D([0], [0], linestyle="-", color=colors[i], label=rf"{co[i]} atm")
            )
        plt.ylabel(r"Current density [Acm$^{-2}$]")
        plt.xlabel("Overpotential [V]")
        plt.yscale("log")
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        solutions_legend = [
            Line2D(
                [0],
                [0],
                color="#1f77b4",
                lw=2,
                linestyle="-",
                label="Moore et. al. 2017",
            ),
            Line2D(
                [0],
                [0],
                color="#ff7f0e",
                lw=2,
                linestyle="-",
                label=r"$\mu$Elektrodica",
            ),
        ]

        first_legend = plt.legend(handles=pressure_legend, loc="lower right")
        plt.gca().add_artist(first_legend)
        plt.legend(handles=solutions_legend, loc="upper left")
        plt.ylim(1e-5, 1e-1)
        plt.xlim(0.18, 0.6)
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

        self.operation.temperature = 50 + 273.15
        self.data.species.c0_reactants = np.array([0.21, 1])
        vdata = 1.18 - p021[:, 1]
        jdata = p021[:, 0]

        F = Fitter(self.data, vdata, jdata)
        self.data.reactions.ga = F.ga_fit
        self.data.species.g_formation_ads = F.gf_fit
        print('\ng_fit:')
        print('ga: ', self.data.reactions.ga)
        print('g_formation: ', self.data.species.g_formation_ads)

        self.nmelek = Calculator(self.data)

        fname = os.path.join(directory, "MooreOxygen_fit.png")
        for label, (mark, linest, data) in po2.items():
            plt.plot(
                1.18 - data[:, 1],
                data[:, 0],
                marker=mark,
                linestyle="",
                color="#1f77b4",
                label="Moore et. al. 2017",
            )

        plt.plot(
            self.operation.potential,
            abs(self.nmelek.results.j),
            linestyle="-",
            color="tab:green",
            label=r"$\mu$Elektrodica Fit",

        )
        plt.ylabel(r"Current density [Acm$^{-2}$]")
        plt.xlabel("Overpotential [V]")
        plt.yscale("log")
        # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')

        plt.legend(loc="upper left")
        plt.ylim(1e-5, 1e-1)
        plt.xlim(0.18, 0.6)
        plt.minorticks_on()
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.show()


if __name__ == "__main__":
    Oxygen()

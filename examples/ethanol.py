"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        examples/ethanol.py

"""
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from melektrodica import Collector, Calculator, Kpynetic, Writer


# for debugging
# from .Tools import showme
# import sys
# sys.exit()


class Ethanol:
    def __init__(self):
        writer = Writer()
        directory = os.path.join("SanchezMonreal2017Ethanol")
        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions
        self.Kpy = Kpynetic(self.data)
        self.melek = Calculator(self.Kpy, 'melek')

        writer.markdown('ethanol', 'theta', self.melek )

        # CSTR mode
        self.data_cstr = copy.deepcopy(self.data)
        self.data_cstr.parameters.cstr = True
        self.data_cstr.reactions.nux = self.data_cstr.reactions.nuc
        self.data_cstr.parameters.Fv = 3e-4
        self.data_cstr.parameters.Ac = 1
        self.Kpy_cstr = Kpynetic(self.data_cstr)
        self.melek_cstr = Calculator(self.Kpy_cstr, 'cstr')



        ch3choh = np.array(
            [
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
                [4.46133e-1, 1.67220e-7],
            ]
        )
        ch3co = np.array(
            [
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
                [4.46618e-1, 6.05464e-1],
            ]
        )
        co = np.array(
            [
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
                [4.49585e-1, 3.71279e-1],
            ]
        )
        ch3 = np.array(
            [
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
                [4.50623e-1, 1.00027e-14],
            ]
        )
        oh = np.array(
            [
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
                [4.46888e-1, 5.28671e-1],
            ]
        )

        species = {
            "CH$_3$CHOH": ("tab:orange", ch3choh, 0),
            "CH$_3$CO": ("tab:blue", ch3co, 1),
            "CO": ("tab:green", co, 2),
            "CH$_3$": ("tab:red", ch3, 3),
            "OH": ("tab:purple", oh, 4),
        }

        fname = os.path.join(directory, "SanchezEthanol_concentrations.png")
        plt.plot(
            self.operation.potential,
            self.melek_cstr.results.c_reactants[:, 0],
            label=r"CH$_3$CH$_2$OH",
        )
        plt.plot(
            self.operation.potential,
            self.melek_cstr.results.c_products[:, 0],
            label=r"CH$_3$CHO",
        )
        plt.xlabel("Overpotential [V]")
        plt.ylabel(r"$c_j$")
        plt.xlim(0, 1)
        # plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()
        plt.legend()
        plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

        fname = os.path.join(directory, "SanchezEthanol_theta.png")
        fname_zoom = os.path.join(directory, "SanchezEthanol_theta_zoom.png")
        for lim in {True, False}:
            for label, (colors, data, idx) in species.items():
                plt.plot(data[:, 0], data[:, 1], linestyle="-", color=colors)
                plt.plot(
                    self.operation.potential,
                    self.melek.results.theta[:, idx],
                    linestyle="--",
                    color=colors,
                )
                plt.plot(
                    self.operation.potential,
                    self.melek_cstr.results.theta[:, idx],
                    linestyle=":",
                    color=colors,
                )

            plt.xlabel("Overpotential [V]")
            plt.ylabel(r"$\theta_j$")
            # plt.grid(visible=True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
            plt.yscale("log")
            plt.minorticks_on()
            plt.tight_layout()

            species_legend = [
                Line2D([0], [0], color=colors, lw=2, label=rf"{specie}*")
                for specie, (colors, _, _) in species.items()
            ]
            solutions_legend = [
                Line2D(
                    [0],
                    [0],
                    color="grey",
                    lw=2,
                    linestyle="-",
                    label="Sanchez-Monreal et. al. 2017",
                ),
                Line2D(
                    [0],
                    [0],
                    color="grey",
                    lw=2,
                    linestyle="--",
                    label=r"$\mu$Elektrodica",
                ),
                Line2D(
                    [0],
                    [0],
                    color="grey",
                    lw=2,
                    linestyle=":",
                    label=r"$\mu$Elektrodica CSTR",
                ),
            ]
            plt.xlim(0, 1)
            plt.ylim(1e-15, 1.5)
            if lim:
                plt.ylim(1e-3, 1.1)
                plt.legend(handles=solutions_legend, loc="lower center")
                plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
                plt.savefig(fname_zoom, dpi=300, bbox_inches="tight", format="png")
            else:
                first_legend = plt.legend(handles=species_legend, loc="lower left")
                plt.gca().add_artist(first_legend)
                # plt.legend(handles=solutions_legend, loc='lower right')
                plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
                plt.savefig(fname, dpi=300, bbox_inches="tight", format="png")
            plt.show()


if __name__ == "__main__":
    Ethanol()

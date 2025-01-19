"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        examples/methanol.py

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from melektrodica import Collector, Calculator


class Methanol:
    def __init__(self):
        directory = os.path.join("ParedesSalazar2023Methanol")

        self.data = Collector(directory)
        self.operation = self.data.parameters
        self.species = self.data.species
        self.reactions = self.data.reactions

        self.melek = Calculator(self.data)

        plt.plot(
            self.operation.potential,
            self.melek.results.theta,
            label=self.data.species.adsorbed,
        )
        plt.legend()
        plt.yscale("log")
        plt.show()

        plt.plot(
            self.operation.potential,
            self.melek.results.fval,
            label=self.data.species.adsorbed,
        )
        plt.legend()
        plt.yscale("log")
        plt.show()


if __name__ == "__main__":
    Methanol()

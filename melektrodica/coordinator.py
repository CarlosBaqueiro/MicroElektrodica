import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import copy

from .kpynetic import Kpynetic
from .tools import Tool


class Coordinator:
    def __init__(self, kpy):
        self.Kpy = copy.deepcopy(kpy)
        self.data = self.Kpy.data

    def plot_rxn_coords_potential(self, source, target, eta, fname):
        self.upsilon = copy.deepcopy(self.data.reactions.nuc)
        self.species = copy.deepcopy(self.data.species.list)
        self.reactions = copy.deepcopy(self.data.reactions.list)
        self.grafo = self.stoichiometric_graphe(self.upsilon, self.species, self.reactions)
        self.paths = self.get_paths(self.grafo, source, target)
        self.pathways = self.get_pathways(self.grafo, self.paths)
        self.sigma = self.get_sigma(self.grafo, self.paths, self.reactions)
        for p, path in enumerate(self.paths):
            figname = f"{fname}_path_{p}.png"
            self.plotter_rxn_coords_potential(p, path, eta, figname)

    def stoichiometric_graphe(self, upsilon, species, reactions):
        """Build a directed graph from a stoichiometric matrix."""
        graphe = nx.MultiDiGraph()
        species = np.array(species)
        for i, row in enumerate(upsilon):
            reaction = reactions[i]
            reactants = species[np.where(row < 0)].tolist()
            products = species[np.where(row > 0)].tolist()
            coeff = [0, 0]
            for reactant in reactants:
                j = np.where(species == reactant)
                coeff[0] = upsilon[i, j].item()
                for product in products:
                    j = np.where(species == product)
                    coeff[1] = upsilon[i, j].item()
                    graphe.add_edge(reactant, product, reaction=reaction, upsilon=coeff)
        return graphe

    def get_paths(self, graphe, source, target):
        paths = list(nx.all_simple_paths(graphe, source, target))
        return paths

    def get_pathways(self, graphe, paths):
        grafo = copy.deepcopy(graphe)
        pathways = []
        for h, nodes in enumerate(paths):
            pathways.append([])
            for n in range(len(nodes) - 1):
                reactant, product = nodes[n], nodes[n + 1]
                edges = grafo[reactant][product]
                key, data = list(edges.items())[-1]
                reaction, upsilon = data['reaction'], data['upsilon']
                pathways[h].append(reaction)
                if key != 0:
                    grafo.remove_edge(reactant, product, key=key)
        del grafo
        return pathways

    def get_sigma(self, graphe, paths, reactions):
        grafo = copy.deepcopy(graphe)
        sigma = np.zeros((len(paths), len(reactions)))
        for h, nodes in enumerate(paths):
            factor = 1
            for n in range(len(nodes) - 1):
                reactant, product = nodes[n], nodes[n + 1]
                edges = grafo[reactant][product]
                key, data = list(edges.items())[-1]
                reaction, upsilon = data['reaction'], data['upsilon']
                i = reactions.index(reaction)
                sigma[h, i] = factor
                factor *= upsilon[1]
                if key != 0:
                    grafo.remove_edge(reactant, product, key=key)
        del grafo
        return sigma

    def get_energies(self, sigma, reactions, eta=0, zero=0):
        dg = self.Kpy.get_argument(eta) * sigma
        ddg = (dg[0, :] - dg[1, :])
        energies = np.zeros(2 * len(reactions) + 1)
        energies[0] = zero
        for i, reaction in enumerate(reactions):
            idx = self.reactions.index(reaction)
            energies[2 * i + 1] = energies[2 * i] + dg[0, idx]
            energies[2 * i + 2] = energies[2 * i] + ddg[idx]
        energies -= energies[zero]
        return energies

    def format_latex_chemical(sef, species):
        chemicals = []
        for chem in species:
            # Añade subíndices: coloca lo que está después de un número como "_" para subíndice
            formatted = ""
            for char in chem:
                if char.isdigit():
                    formatted += f"_{char}"  # Usa subíndice en LaTeX
                elif char == '+':
                    formatted += "^+"
                else:
                    formatted += char
            chemicals.append(f"$\\mathrm{{{formatted}}}$")  # Agregar delimitadores de LaTeX
        return chemicals

    def plotter_rxn_coords_potential(self, p, path, potential, figname):
        fig, ax = plt.subplots()
        plateau = 1 / 4
        reactions = self.pathways[p]
        species = Tool.format_latex_chemical(path)
        potential_legend = []
        colors = plt.cm.viridis(np.linspace(0, 1, len(potential)))
        colors[0] = plt.cm.tab10(1)
        for j, eta in enumerate(potential):
            energies = self.get_energies(self.sigma[p], reactions, eta)
            print(f'eta = {eta} V, energies = {energies} eV,')
            color = colors[j]
            potential_legend.append(
                Line2D([0], [0], linestyle="-", color=color, label=rf"$\eta = {eta}$ V")  # Assuming eta is a potential
            )
            for i, g in enumerate(energies):
                plt.plot([(i + 1 - plateau), (i + 1 + plateau)], [g, g], color=color, linewidth=2)
                if i <= len(energies) - 2:
                    plt.plot([(i + 1 + plateau), (i + 2 - plateau)], [g, energies[i + 1]],
                             linestyle=":", color=color, linewidth=0.5)
        x_label = [""] * len(energies)
        x_label[0::2] = species
        x_label[1::2] = reactions
        plt.plot([1 + plateau, 2 * len(reactions) + 1 + plateau], [0, 0], color='gray', linestyle="--", linewidth=1)
        plt.xlim(1 - plateau, 2 * len(reactions) + 1 + plateau)
        plt.xlabel("Reaction Coordinate")
        plt.ylabel(r"$\Delta G_{r,\xi_h}$ [eV]")
        plt.xticks(range(1, len(energies) + 1), x_label)
        plt.legend(handles=potential_legend)
        plt.minorticks_on()
        #ax.tick_params(axis='both', which='both', direction='in')
        ax.tick_params(axis='x', which='minor', bottom=False)
        plt.tight_layout()
        plt.savefig(figname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

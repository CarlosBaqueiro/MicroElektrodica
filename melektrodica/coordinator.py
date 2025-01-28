#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    μElektrodica© 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Coordinator class

"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import copy

from .kpynetic import Kpynetic
from .tools import Tool


class Coordinator:
    """
    Represents a coordinator for analyzing reaction coordinates and potentials within
    a chemical or biochemical network.

    Provides functionality to construct stoichiometric graphs, determine pathways between
    species, compute energy changes along reaction pathways, and visualize reaction coordinates
    and potential energy surfaces.

    Attributes
    ----------
    Kpy : copy.deepcopy
        Deep copy of the input data object containing information about reactions, species,
        and other necessary attributes for energy and pathway computations.
    data : object
        Direct access to the data contained within `Kpy`.
    """

    def __init__(self, kpy):
        """
        Instantiates a new object using a deep copy of the provided `kpy` object. The copied
        object's data attribute is assigned directly to the instance's `data` attribute. This
        ensures the immutability of the provided input and initializes the data attribute
        for further operations.

        Parameters
        ----------
        kpy : object
            The object to be deeply copied and used for initialization. Its `data` attribute
            will be accessed and assigned to the instance.
        """
        self.Kpy = copy.deepcopy(kpy)
        self.data = self.Kpy.data

    def plot_rxn_coords_potential(self, source, target, eta, fname):
        """
        Plots reaction coordinates and potential energy for pathways between the source and the target.

        This function primarily utilizes information from the stoichiometric graph to identify
        reaction pathways and compute reaction coordinates and potentials. It produces individual
        plots for each identified pathway, illustrating the progress of reactions and corresponding
        potential energy changes.

        Parameters
        ----------
        source : str
            Identifier for the starting species in the pathway.
        target : str
            Identifier for the target species in the pathway.
        eta : float
            Parameter affecting the visualization of reaction coordinates and potential energy. Its exact
            purpose depends on the implementation of the plotting function.
        fname : str
            Base file name to save the generated plots. Each pathway will be saved as an individual
            file, with a unique suffix indicating the pathway index.
        """
        self.upsilon = copy.deepcopy(self.data.reactions.upsilon_c)
        self.species = copy.deepcopy(self.data.species.list)
        self.reactions = copy.deepcopy(self.data.reactions.list)
        self.grafo = self.stoichiometric_graphe(
            self.upsilon, self.species, self.reactions
        )
        self.paths = self.get_paths(self.grafo, source, target)
        self.pathways = self.get_pathways(self.grafo, self.paths)
        self.sigma = self.get_sigma(self.grafo, self.paths, self.reactions)
        for p, path in enumerate(self.paths):
            figname = f"{fname}_path_{p}.png"
            self.plotter_rxn_coords_potential(p, path, eta, figname)

    def stoichiometric_graphe(self, upsilon, species, reactions):
        """
        Constructs a stoichiometric graph that represents the transformation of species through
        a series of reactions. The graph edges represent the mapping of reactants to products
        linked by a specific reaction and associated stoichiometric coefficients.

        Parameters:
            upsilon: numpy.ndarray
                A 2D array representing stoichiometric coefficients for the reactions
                with rows corresponding to reactions and columns to species.
            species: list or array-like
                An ordered collection of species identifiers corresponding to the
                columns of the upsilon matrix.
            reactions: list
                A collection of identifiers for reactions corresponding to the rows
                of the upsilon matrix.

        Returns:
            networkx.MultiDiGraph
                A directed multigraph where nodes represent species and edges represent
                the reactions. Edges are annotated with reaction IDs and their respective
                stoichiometric coefficients.

        Raises:
            None
        """
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
        """
        Find all simple paths in a graph from the source node to the target node.

        This function utilizes the NetworkX library to compute all simple paths
        between a source node and a target node in a given graph. Simple paths are
        paths that do not revisit any node within the graph, ensuring no cycles are present.

        Parameters
        ----------
        graphe : networkx.Graph
            The input graph for which paths need to be calculated. This can be any
            graph object supported by NetworkX, such as Graph, DiGraph, MultiGraph,
            or MultiDiGraph.
        source : Hashable
            The starting node of the paths. This must be a node present within the
            graph.
        target : Hashable
            The ending node of the paths. This must be a node present within the
            graph.

        Returns
        -------
        list of list
            A list of paths, where each path is represented as a list of nodes that
            constitute the path from the source to the target node.
        """
        paths = list(nx.all_simple_paths(graphe, source, target))
        # Filter to prevent "short circuits"
        for h, nodes in enumerate(paths):
            for n in range(len(nodes) - 1):
                reactant, product = nodes[n], nodes[n + 1]
                if reactant == source and product == target:
                    graphe.remove_edge(reactant, product)
        paths = list(nx.all_simple_paths(graphe, source, target))
        return paths

    def get_pathways(self, graphe, paths):
        """
        Generate pathways based on the provided graph and paths.

        This function creates deep copies of the input graph data to avoid modifying
        the original graph. It processes the paths to extract the sequence of reactions
        for each path and removes certain edges based on their key values. Finally,
        it returns the list of pathways comprising the ordered reactions.

        Parameters
        ----------
        graphe : networkx.Graph
            A directed graph representing the system, where nodes represent entities
            and edges represent reactions between them. Edge data must include a
            'reaction' (str or identifier) and 'upsilon' (any relevant value) field.
        paths : list of list
            A list of paths, where each path is represented as a list of nodes (vertex
            identifiers) forming a valid sequence in the given graph.

        Returns
        -------
        list of list
            A 2D list, where the outer list represents the collection of paths, and
            each inner list contains the ordered reactions corresponding to a single
            path in the input graph.

        """
        grafo = copy.deepcopy(graphe)
        pathways = []
        for h, nodes in enumerate(paths):
            pathways.append([])
            for n in range(len(nodes) - 1):
                reactant, product = nodes[n], nodes[n + 1]
                edges = grafo[reactant][product]
                key, data = list(edges.items())[-1]
                reaction, upsilon = data["reaction"], data["upsilon"]
                pathways[h].append(reaction)
                if key != 0:
                    grafo.remove_edge(reactant, product, key=key)
        del grafo
        return pathways

    def get_sigma(self, graphe, paths, reactions):
        """
        Computes the sigma matrix for the given graph, paths, and reactions.

        This function calculates the sigma matrix, which represents the
        contributions of each reaction in the given paths, using a deepcopy of the
        provided graph. It iterates through the paths, processes each pair of
        nodes, and updates the sigma matrix accordingly by considering reaction
        and upsilon data. Additionally, edges in the graph are removed based on
        specific keys during the computation.

        Parameters
        ----------
        graphe : dict
            A graph represented as a dictionary where nodes and edges with
            associated data (reaction and upsilon) are specified.
        paths : list of list of any
            A list containing the paths, where each path is a list of nodes
            (reactants and products) in order.
        reactions : list of str
            A list of all the reactions to be considered, where each reaction is
            represented as a string.

        Returns
        -------
        numpy.ndarray
            A 2D array (sigma matrix) where each row corresponds to a path and
            each column corresponds to a reaction. The values represent the
            contribution of each reaction within a path.
        """
        grafo = copy.deepcopy(graphe)
        sigma = np.zeros((len(paths), len(reactions)))
        for h, nodes in enumerate(paths):
            factor = 1
            for n in range(len(nodes) - 1):
                reactant, product = nodes[n], nodes[n + 1]
                edges = grafo[reactant][product]
                key, data = list(edges.items())[-1]
                reaction, upsilon = data["reaction"], data["upsilon"]
                i = reactions.index(reaction)
                sigma[h, i] = factor
                factor *= upsilon[1]
                if key != 0:
                    grafo.remove_edge(reactant, product, key=key)
        del grafo
        return sigma

    def get_energies(self, sigma, reactions, eta=0, zero=0):
        """
        Calculate the energy profile for given reactions.

        This method computes the energy profile of chemical reactions based on the
        input parameters. It utilizes the free energy changes associated with the
        reaction pathways to produce an array containing energy levels for each
        reaction step including intermediates.

        Parameters
        ----------
        sigma : ndarray
            Free energy changes, typically as a 2D array where rows correspond to
            different reaction pathways.
        reactions : list of str
            List of reaction identifiers corresponding to the reactions for which
            the energy profile is to be calculated.
        eta : int, optional
            A parameter used to compute the free energy changes (default is 0). It
            affects the argument passed to the internal function `get_argument`.
        zero : int, optional
            The index of the energy level used as the reference energy (default is 0).

        Returns
        -------
        energies : ndarray
            A 1D array containing the energy levels for each reaction step and
            intermediate states. The energy levels are shifted such that the
            `zero`-indexed level is at 0.

        """
        dg = self.Kpy.get_argument(eta) * sigma
        ddg = dg[0, :] - dg[1, :]
        energies = np.zeros(2 * len(reactions) + 1)
        energies[0] = zero
        for i, reaction in enumerate(reactions):
            idx = self.reactions.index(reaction)
            energies[2 * i + 1] = energies[2 * i] + dg[0, idx]
            energies[2 * i + 2] = energies[2 * i] + ddg[idx]
        energies -= energies[zero]
        return energies

    def plotter_rxn_coords_potential(self, p, path, potential, figname):
        """
        Plots reaction coordinates against potential energies for specified pathway and potential values.

        This method visualizes the energetic profile of a chemical reaction pathway as a
        function of reaction coordinates (`path`) under varying applied potentials (`potential`).
        The produced plot assists in understanding the free energy landscape for reactions under
        different electrochemical conditions, with reaction steps and intermediates represented
        along the x-axis, and their corresponding free energy changes shown on the y-axis.
        Potential-dependent energetic variations are color-coded and distinguished by unique
        legends.

        Parameters
        ----------
        p : int
            Index specifying the pathway to be analyzed within the object's pathways.

        path : list
            List of chemical species or reaction intermediates involved in the reaction pathway.

        potential : list of float
            The set of applied potentials at which the energy profiles are evaluated.

        figname : str
            The name of the file where the generated plot is to be saved.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If inputs are not properly formatted or incompatible with the plot generation.
        FileNotFoundError
            If the specified file path for saving the plot is invalid.
        """
        fig, ax = plt.subplots()
        plateau = 1 / 4
        reactions = self.pathways[p]
        species = Tool.format_latex_chemical(path)
        potential_legend = []
        colors = plt.cm.viridis(np.linspace(0, 1, len(potential)))
        colors[0] = plt.cm.tab10(1)
        for j, eta in enumerate(potential):
            energies = self.get_energies(self.sigma[p], reactions, eta)
            # print(f'eta = {eta} V, energies = {energies} eV,')
            color = colors[j]
            potential_legend.append(
                Line2D(
                    [0], [0], linestyle="-", color=color, label=rf"$\eta = {eta}$ V"
                )  # Assuming eta is a potential
            )
            for i, g in enumerate(energies):
                plt.plot(
                    [(i + 1 - plateau), (i + 1 + plateau)],
                    [g, g],
                    color=color,
                    linewidth=2,
                )
                if i <= len(energies) - 2:
                    plt.plot(
                        [(i + 1 + plateau), (i + 2 - plateau)],
                        [g, energies[i + 1]],
                        linestyle=":",
                        color=color,
                        linewidth=0.5,
                    )
        x_label = [""] * len(energies)
        x_label[0::2] = species
        x_label[1::2] = reactions
        plt.plot(
            [1 + plateau, 2 * len(reactions) + 1 + plateau],
            [0, 0],
            color="gray",
            linestyle="--",
            linewidth=1,
        )
        plt.xlim(1 - plateau, 2 * len(reactions) + 1 + plateau)
        plt.xlabel("Reaction Coordinate")
        plt.ylabel(r"$\Delta G_{r,\xi_h}$ [eV]")
        plt.xticks(range(1, len(energies) + 1), x_label)
        plt.legend(handles=potential_legend)
        plt.minorticks_on()
        # ax.tick_params(axis='both', which='both', direction='in')
        ax.tick_params(axis="x", which="minor", bottom=False)
        plt.tight_layout()
        plt.savefig(figname, dpi=300, bbox_inches="tight", format="png")
        plt.show()

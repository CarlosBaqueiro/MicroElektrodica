"""

    μElektrodica (Uxmal, version 1.0.0)
        A Python Tool for Modeling Microkinetic Electrocatalytic Reactions
        Copyright (C) 2024 C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        All rights reserved.

        Collector class

"""

import os
import numpy as np
import re

from .Tools import unit_conversion

# for debugging
#import sys
# sys.exit()
#from .Tools import showme

class DataParameters:
    """
    A class used to manage and initialize simulation parameters for different chemical and
    experimental scenarios. It processes the parameters and variables from the input file,
    performs unit conversions, and initializes various system components like anode, CSTR,
    and reaction parameters.

    :ivar parameters_list: A list containing the parameter names extracted from the file.
    :type parameters_list: numpy.ndarray
    :ivar variables_list: A list containing the variable names extracted from the file.
    :type variables_list: numpy.ndarray
    :ivar potential: A numpy array representing the range of potential values.
    :type potential: numpy.ndarray
    :ivar T: Temperature in Kelvin after unit conversion.
    :type T: float
    :ivar anode: Electrode, True for anode, False for cathode.
    :type anode: bool
    :ivar cstr: Initialized Continuous Stirred-Tank Reactor (CSTR) parameters.
    :type cstr: bool
    :ivar tst: Initialized Transient State Theory parameters.
    :type tst: bool
    :ivar js: Initialized j* parameters.
    :type js: bool
    :ivar experimental: Initialized experimental parameters.
    :type experimental: bool
    :ivar chemical: Initialized chemical parameters.
    :type chemical: bool
    :ivar Fv: Volumetric flux value for CSTR, if applicable.
    :type Fv: float
    :ivar Ac: Catalyst active surface area for CSTR, if applicable.
    :type Ac: float
    :ivar pre_exponential: Pre-exponential factor in rate constant calculation.
    :type pre_exponential: float
    :ivar js_value: Value of j*, if applicable.
    :type js_value: float
    :ivar kappa: Kappa value for transient state theory, if applicable.
    :type kappa: float
    :ivar m: m value for transient state theory, if applicable.
    :type m: float
    :ivar DG_reaction: DG reaction parameter for chemical data, if applicable.
    :type DG_reaction: bool
    :ivar G_formation: G formation parameter for chemical data, if applicable.
    :type G_formation: bool
    """

    def __init__(self, parameters_file: str) -> None:

        """
            A class used to manage and initialize simulation parameters for different chemical
            and experimental scenarios. It processes the parameters and variables from the
            input file, performs unit conversions, and initializes various system components
            like anode, CSTR, and reaction parameters.

            :param parameters_file: The file containing the simulation parameters and variables.
            :type parameters_file: str
            """
        header, raw_data = Collector.raw_data(parameters_file)
        self.parameters_list = raw_data[:, header.index('Parameters')]
        self.variables_list = raw_data[:, header.index('Variables')]
        values = raw_data[:, header.index('Value')]
        units = raw_data[:, header.index('Units')]

        # Temperature
        T = float(values[self.variables_list == 'Temperature'])
        uT = units[self.variables_list == 'Temperature']
        self.T = unit_conversion('Temperature', T, uT, 'K')

        # Potential
        iE = float(values[self.variables_list == 'Initial potential'])
        fE = float(values[self.variables_list == 'Final potential'])
        hE = float(values[self.variables_list == 'Step potential'])
        self.potential = np.arange(iE, fE + hE, hE) 


        # Initialize parameters
        self.anode = initialize(values, self.variables_list, 'Anode')
        self.cstr = initialize(values, self.parameters_list, 'Continuous Stirred-Tank Reactor')
        self.tst = initialize(values, self.parameters_list, 'Transient state theory')
        self.js = initialize(values, self.variables_list, 'j*')
        self.experimental = initialize(values, self.parameters_list, 'Experimental')
        self.chemical = initialize(values, self.parameters_list, 'Chemical')
        

        # Continuous Stirred-Tank Reactor (CSTR): Concentrations in function of potential
        if self.cstr:
            self.Fv = float(values[self.variables_list == 'Volumetric flux'])
            self.Ac = float(values[self.variables_list == 'Catalyst Active surface area'])

        # Rate Constants
        # Pre-exponential
        self.pre_exponential = float(values[self.variables_list == 'A'])
        if self.js:
            self.js_value = float(values[self.variables_list == 'j* (value)'])
        if self.tst:
            self.kappa = float(values[self.variables_list == 'kappa'])
            self.m = float(values[self.variables_list == 'm'])

        # Check Chemical part details
        if self.chemical:
            self.DG_reaction = initialize(values, self.variables_list, 'DG_reaction')
            self.G_formation = initialize(values, self.variables_list, 'G_formation')

    # TODO: Add data recollection for more models and operations conditions

def initialize(values: list, lista: list, name: str) -> bool:
    """
    Checks if a given name exists in a list and if the corresponding value in another
    list is 'True'.

    :param values: List of values where each position corresponds to an element
                   in 'lista'.
    :param lista: List of names to be checked.
    :param name: Name to be searched within 'lista'.

    :return: Returns True if 'name' is found in 'lista' and the corresponding value
             in 'values' is 'True'. Otherwise, returns False.
    """
    if name in lista and values[lista == name] == 'True':
        return True
    else:
        return False

class DataSpecies:
    """
    Represents data related to chemical species, including reactants,
    products, adsorbed species, and catalysts.

    Detailed description of the DataSpecies class that initializes and
    processes information related to various chemical species from a given
    data file. It categorizes species based on their roles (such as reactants,
    products, adsorbed, and catalysts) and computes relevant properties,
    including formation energies and initial concentrations.

    :ivar reactants: List of reactant species.
    :type reactants: list of str
    :ivar products: List of product species.
    :type products: list of str
    :ivar adsorbed: List of adsorbed species.
    :type adsorbed: list of str
    :ivar catalyst: List of catalyst species.
    :type catalyst: list of str
    :ivar c0_reactants: Initial concentrations of reactant species.
    :type c0_reactants: numpy.ndarray
    :ivar c0_products: Initial concentrations of product species.
    :type c0_products: numpy.ndarray
    :ivar list: Combined list of all species including an electron placeholder.
    :type list: list of str
    :ivar nu_catalyst: Matrix indicating presence of adsorbed species in each catalyst.
    :type nu_catalyst: numpy.ndarray
    :ivar G_formation_rct: Gibbs free energy of formation for reactants.
    :type G_formation_rct: numpy.ndarray
    :ivar G_formation_ads: Gibbs free energy of formation for adsorbed species.
    :type G_formation_ads: numpy.ndarray
    :ivar G_formation_prd: Gibbs free energy of formation for products.
    :type G_formation_prd: numpy.ndarray
    """
    def __init__(self,species_file, parameters):
        header, raw_data = Collector.raw_data(species_file)
        species_list = raw_data[:,header.index('Species')]
        index = header.index('RPACe')
        self.reactants = species_list[raw_data[:, index] == 'R'].tolist()
        self.products  = species_list[raw_data[:, index] == 'P'].tolist()
        self.adsorbed = species_list[raw_data[:, index] == 'A'].tolist()
        self.catalyst   = species_list[raw_data[:, index] == 'C'].tolist()
        self.c0_reactants = np.array(raw_data[:, header.index('c0')][raw_data[:, index] == 'R'].astype(float))
        self.c0_products  = np.array(raw_data[:, header.index('c0')][raw_data[:, index] == 'P'].astype(float))
        self.list = self.reactants + self.products + self.adsorbed + self.catalyst + ['e-']

        self.nu_catalyst = np.zeros((len(self.catalyst), len(self.adsorbed)))
        for i in range(len(self.catalyst)):
            speciesincatalyst = species_list[raw_data[:, header.index('Catalyst')] == self.catalyst[i]]
            for specie in speciesincatalyst:
                if specie in self.adsorbed:
                    self.nu_catalyst[i, self.adsorbed.index(specie)] = 1

        if parameters.chemical:
            if parameters.G_formation:
                self.G_formation_rct = np.array(raw_data[:,header.index('DG_formation')][raw_data[:, index] == 'R'].astype(float))
                self.G_formation_ads = np.array(raw_data[:,header.index('DG_formation')][raw_data[:, index] == 'A'].astype(float))
                self.G_formation_prd = np.array(raw_data[:,header.index('DG_formation')][raw_data[:, index] == 'P'].astype(float))


class DataReactions:
    """
    Manages and processes reaction data.

    This class reads reaction data from a file, processes it to extract various
    parameters and coefficients, and makes the data accessible for further
    analysis. It also handles different experimental and chemical conditions
    based on the provided parameters.

    :ivar list: List of reaction IDs.
    :type list: list
    :ivar beta: Array of Beta parameters for the reactions.
    :type beta: numpy.ndarray
    :ivar nu: Stoichiometric matrix for the reactions.
    :type nu: numpy.ndarray
    :ivar ne: Array of the number of electrons transferred in the reactions.
    :type ne: numpy.ndarray
    :ivar nuc: Stoichiometric matrix without catalysts.
    :type nuc: numpy.ndarray
    :ivar nua: Stoichiometric matrix for adsorbates.
    :type nua: numpy.ndarray
    :ivar nux: Stoichiometric matrix adjusted for CSTR or experimental conditions.
    :type nux: numpy.ndarray
    :ivar k_f: Array of forward reaction rate constants (optional).
    :type k_f: numpy.ndarray
    :ivar k_b: Array of backward reaction rate constants (optional).
    :type k_b: numpy.ndarray
    :ivar Ga: Array of Gibbs free energy changes (optional).
    :type Ga: numpy.ndarray
    :ivar DG_reaction: Array of Delta G reaction values (optional).
    :type DG_reaction: numpy.ndarray
    """
    def __init__(self, reaction_file, parameters, species):
        header, raw_data = Collector.raw_data(reaction_file)
        self.list = raw_data[:, header.index('id')].tolist()
        self.beta = np.array(raw_data[:, header.index('Beta')].astype(float))

        self.nu = np.zeros((len(self.list), len(species.list)))
        for i in range(len(self.list)):
            r = raw_data[:, header.index('Reactions')][i]
            left, right = re.split(r'<->', r)

            speciesinreaction, stoichiometric = self.process_reaction(left, species.list)
            for specie, coeff in zip(speciesinreaction, stoichiometric):
                self.nu[i, species.list.index(specie)] = -coeff

            speciesinreaction, stoichiometric = self.process_reaction(right, species.list)
            for specie, coeff in zip(speciesinreaction, stoichiometric):
                self.nu[i, species.list.index(specie)] = coeff

        self.ne = self.nu[:,-1]                         # Number of electrons transferred
        self.nu = self.nu[:,:-1]                        # All coefficients, catalysts included
        self.nuc = self.nu[:,:-len(species.catalyst)]   # All coefficients, without catalysts
        self.nua = self.nuc[:,-len(species.adsorbed):]    # Adsorbates coefficients

        if parameters.cstr:
            self.nux = self.nuc
        else:
            self.nux = self.nua

        if parameters.experimental:
            self.k_f = np.array(raw_data[:, header.index('k_f')].astype(float))
            self.k_b = np.array(raw_data[:, header.index('k_b')].astype(float))

        if parameters.chemical:
            self.Ga = np.array(raw_data[:, header.index('Ga')].astype(float))
            if parameters.DG_reaction:
                self.DG_reaction = np.array(raw_data[:, header.index('DG_reaction')].astype(float))

        #TODO: Add data recollection for more models

    @staticmethod
    def process_reaction(side, species_list):
        speciesinreaction =[]
        stoichiometric = []
        str = re.split(r'\s+', side)
        for s in str:
            s = s.strip()
            if s in ['', '+']:
                continue
            match = re.search(r"(\d+(\.\d+)?)?(.+)", s)
            if match:
                coefficient = float(match.group(1)) if match.group(1) else 1
                specie = match.group(3)
                if specie in species_list:
                    speciesinreaction.append(specie)
                    stoichiometric.append(float(coefficient))
                else:
                    raise ValueError(f"Process reaction error:  {specie} not found")
            else:
                raise ValueError(f"Process reaction error:  {str} not found")
        return speciesinreaction, stoichiometric


class Collector:
    """
    Manages and processes data related to parameters, species, and reactions from specified directory.

    This class is designed to handle and manage data from files located in a given directory.
    It initializes various data parameters, species, and reactions by reading from corresponding files.

    :ivar parameters: Instance handling data parameters from the parameters.md file.
    :type parameters: DataParameters
    :ivar species: Instance managing species information from the species.md file.
    :type species: DataSpecies
    :ivar reactions: Instance managing reaction data from the reactions.md file.
    :type reactions: DataReactions
    """
    def __init__(self, directory):
        self.parameters = DataParameters(os.path.join(directory, 'parameters.md'))
        self.species = DataSpecies(os.path.join(directory, 'species.md'), self.parameters)
        self.reactions = DataReactions(os.path.join(directory, 'reactions.md'), self.parameters, self.species)

    @staticmethod
    def raw_data(name_file):
        """
            Processes raw data from a file with a specific format.

            The class contains methods to read and process data from a specified file.
            """
        if not os.path.exists(name_file):
            raise FileNotFoundError(f"File was not found at the expected location: {name_file}")
        print(f'\t{name_file}')
        header = []
        raw_data = []
        with open(name_file, 'r') as f:
            lines = f.readlines()
            header = [col.strip() for col in lines[0].split('|')[1:-1]]
            for line in lines[2:]:
                if line.strip() and line[0] != '#':
                    entry_data = [col.strip() for col in line.split('|')[1:-1]]
                    raw_data.append(entry_data)
        raw_data = np.array(raw_data)
        return header, raw_data




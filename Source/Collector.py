"""
     μElektrodica:
                     Collector class documentation
     @Author : Carlos Baqueiro Basto

"""
import os
import numpy as np
import re

from .Tools import unit_conversion

class DataParameters:
    """
    A class used to manage and initialize simulation parameters for different chemical
    and experimental scenarios. It processes the parameters and variables from the
    input file, performs unit conversions, and initializes various system components
    like anode, CSTR, and reaction parameters.

    :param parameters_file: The file containing the simulation parameters and variables.
    :type parameters_file: str

    :ivar parameters_list: List of parameter names extracted from the input file.
    :ivar variables_list: List of variable names extracted from the input file.
    :ivar T: Temperature in Kelvin, after unit conversion.
    :ivar potential: Array of potential values generated from initial, final, and step potential.
    :ivar anode: Initialization status of the anode component.
    :ivar cstr: Initialization status of the Continuous Stirred-Tank Reactor component.
    :ivar tst: Initialization status of the Transient State Theory component.
    :ivar experimental: Initialization status of the experimental component.
    :ivar chemical: Initialization status of the chemical component.
    :ivar Fv: Volumetric flux, initialized if CSTR component is true.
    :ivar Ac: Catalyst active surface area, initialized if CSTR component is true.
    :ivar kappa: Pre-exponential rate constant kappa, initialized if TST component is true.
    :ivar m: Pre-exponential rate constant m, initialized if TST component is true.
    :ivar pre_exponential: Pre-exponential rate constant, initialized if TST component is false.
    :ivar DG_reaction: Gibbs free energy of reaction, initialized if chemical component is true.
    :ivar G_formation: Gibbs free energy of formation, initialized if chemical component is true.
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
        self.parameters_list: list = raw_data[:, header.index('Parameters')]
        self.variables_list: list = raw_data[:, header.index('Variables')]
        values: list = raw_data[:, header.index('Value')]
        units: list = raw_data[:, header.index('Units')]

        # Temperature
        T: float = float(values[self.variables_list == 'Temperature'])
        uT: str = units[self.variables_list == 'Temperature']
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
        self.experimental = initialize(values, self.parameters_list, 'Experimental')
        self.chemical = initialize(values, self.parameters_list, 'Chemical')
        

        # Continuous Stirred-Tank Reactor (CSTR): Concentrations in function of potential
        if self.cstr:
            self.Fv = float(values[self.variables_list == 'Volumetric flux'])
            self.Ac = float(values[self.variables_list == 'Catalyst Active surface area'])

        # Rate Constants
        # Pre-exponential
        if self.tst:
            self.kappa = float(values[self.variables_list == 'kappa'])
            self.m = float(values[self.variables_list == 'm'])
        else:
            self.pre_exponential = float(values[self.variables_list == 'j* (or A)'])

        # Check Chemical part details
        if self.chemical:
            self.DG_reaction = initialize(values, self.variables_list, 'DG_reaction')
            self.G_formation = initialize(values, self.variables_list, 'G_formation')

    # TODO: Add data recollection for more models and operations conditions

def initialize(values: list, lista: list, name: str) -> bool:
    """
    Initialize the system based on the provided values, list, and name. This function checks if the name is present in
    the list and if the corresponding value in the list is 'True'. If both conditions are met, the function returns True.
    Otherwise, it returns False.

    :param lista:
    :param values: Dictionary with boolean-like values mapped to lists.
    :type values: dict
    :param list: List of names to be checked.
    :type list: list
    :param name: Name to be checked in the list.
    :type name: str
    :return: True if the name is in the list and the corresponding value is 'True', otherwise False.
    :rtype: bool
    """
    if name in lista and values[lista == name] == 'True':
        return True
    else:
        return False

class DataSpecies:
    """
    Handles species data, including reactants, products, adsorbed species, and catalysts,
    as well as their initial concentrations and formation Gibbs energies. It utilizes raw
    data to categorize species and compile comprehensive lists and matrices useful for
    further chemical computations and simulations.

    :ivar reactants: List of reactant species.
    :type reactants: list
    :ivar products:  list of product species.
    :type products: list
    :ivar adsorbed: list of adsorbed species.
    :type adsorbed: list
    :ivar catalyst: list  of catalyst species.
    :type catalyst: list
    :ivar c0_reactants: Initial concentrations of reactant species.
    :type c0_reactants: numpy.ndarray
    :ivar c0_products: Initial concentrations of product species.
    :type c0_products: numpy.ndarray
    :ivar list: Comprehensive list of all species including an electron.
    :type list: list
    :ivar nu_catalyst: Matrix indicating which adsorbed species are in which catalyst.
    :type nu_catalyst: numpy.ndarray
    :ivar G_formation_rct: Formation Gibbs energies of reactant species.
    :type G_formation_rct: numpy.ndarray
    :ivar G_formation_ads: Formation Gibbs energies of adsorbed species.
    :type G_formation_ads: numpy.ndarray
    :ivar G_formation_prd: Formation Gibbs energies of product species.
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
    This class deals with biochemical reactions data processing.

    Detailed processing for reactions, including splitting into reactants and products,
    calculating stoichiometric coefficients, and differentiating between various reaction
    parameters based on provided configurations.

    :ivar list: List of reaction ids.
    :type list: list
    :ivar beta: Beta coefficients for reactions.
    :type beta: numpy.ndarray
    :ivar nu: Stoichiometric coefficients matrix.
    :type nu: numpy.ndarray
    :ivar ne: Number of electrons transferred in reactions.
    :type ne: numpy.ndarray
    :ivar nuc: Stoichiometric coefficients excluding catalysts.
    :type nuc: numpy.ndarray
    :ivar nua: Coefficients specific to adsorbates.
    :type nua: numpy.ndarray
    :ivar nux: Coefficients based on whether CSTR or not.
    :type nux: numpy.ndarray
    :ivar k_f: Forward rate constants (if experimental).
    :type k_f: numpy.ndarray
    :ivar k_b: Backward rate constants (if experimental).
    :type k_b: numpy.ndarray
    :ivar Ga: Gibbs free energy changes (if chemical).
    :type Ga: numpy.ndarray
    :ivar DG_reaction: Change in free energy for reaction (if provided).
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
    Handles the collection and processing of data from a specified directory.

    The class initializes with the directory path and reads parameters, species,
    and reactions data from associated files.

    :ivar parameters: Contains parameters data read from 'parameters.md'.
    :type parameters: DataParameters
    :ivar species: Contains species data read from 'species.md'.
    :type species: DataSpecies
    :ivar reactions: Contains reactions data read from 'reactions.md'.
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




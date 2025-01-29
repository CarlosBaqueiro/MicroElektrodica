"""

    μElektrodica © 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Collector class

"""

import os
import numpy as np
import re
from .writer import Writer
from .tools import Tool


# for debugging
# import sys
# sys.exit()


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
    :ivar temperature: Temperature in Kelvin after unit conversion.
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

    def __init__(self, parameters_file: str, writer: object) -> None:
        """
        A class used to manage and initialize simulation parameters for different chemical
        and experimental scenarios. It processes the parameters and variables from the
        input file, performs unit conversions, and initializes various system components
        like anode, CSTR, and reaction parameters.

        :param parameters_file: The file containing the simulation parameters and variables.
        :type parameters_file: str
        """
        writer.message(f"Reading parameters from file: {parameters_file}")
        header, raw_data = Collector.raw_data(parameters_file, writer=writer)
        self.parameters_list = raw_data[:, header.index("Parameters")]
        self.variables_list = raw_data[:, header.index("Variables")]
        values = raw_data[:, header.index("Value")]
        units = raw_data[:, header.index("Units")]
        writer.message(f"Parameters:")
        self.anode = initialize(values, self.variables_list, "Anode")
        writer.message(f"Anode: {self.anode}")

        # Temperature
        temparature = float(values[self.variables_list == "Temperature"])
        temparature_units = units[self.variables_list == "Temperature"]
        self.temperature = Tool.unit_conversion(
            "Temperature", temparature, temparature_units, "K"
        )

        # Potential
        _initial = float(values[self.variables_list == "Initial potential"])
        _final = float(values[self.variables_list == "Final potential"])
        _increment = float(values[self.variables_list == "Step potential"])
        self.potential = np.arange(_initial, _final + _increment, _increment)

        writer.message(
            f"Operation conditions:\n"
            f"\t\tTemperature: {self.temperature} K\n"
            f"\t\tOverpotential: [{_initial}: {_increment}: {_final}] V"
        )

        # Initialize parameters

        self.tst = initialize(values, self.parameters_list, "Transient state theory")
        self.js = initialize(values, self.variables_list, "j*")
        self.experimental = initialize(values, self.parameters_list, "Experimental")
        self.thermochemical = initialize(values, self.parameters_list, "Thermochemical")
        # Rate Constants
        # Pre-exponential
        self.pre_exponential = float(values[self.variables_list == "A"])
        message = "Pre-exponential factor: "
        if self.js:
            self.js_value = float(values[self.variables_list == "j* (value)"])
            message += f"{self.js_value} (from j*)"
        elif self.tst:
            self.kappa = float(values[self.variables_list == "kappa"])
            self.m = float(values[self.variables_list == "m"])
            message += f"kappa: {self.kappa}, m: {self.m} (from TST)"
        else:
            message += f"{self.pre_exponential} (default value)"
        writer.message(message)

        if self.experimental:
            writer.message(f"Experimental kinetics rate constants: {self.experimental}")

        # Check Chemical part details
        if self.thermochemical:
            message = "Thermochemical part details:\n"
            self.dg_reaction = initialize(values, self.variables_list, "DG_reaction")
            self.g_formation = initialize(values, self.variables_list, "G_formation")
            message += f"\tDG_reaction: {self.dg_reaction}\n"
            message += f"\tG_formation: {self.g_formation}"
            writer.message(message)
        # Continuous Stirred-Tank Reactor (CSTR): Concentrations in function of potential
        self.cstr = initialize(
            values, self.parameters_list, "Continuous Stirred-Tank Reactor"
        )
        if self.cstr:
            self.volumetric_flow = float(
                values[self.variables_list == "Volumetric flux"]
            )
            self.catalyst_area = float(
                values[self.variables_list == "Catalyst Active surface area"]
            )
            writer.message(
                f"Continuous Stirred-Tank Reactor model: {self.cstr}\n"
                f"\tVolumetric flux: {self.volumetric_flow}\n"
                f"\tCatalyst Active surface area: {self.catalyst_area}"
            )

    # TODO: Add data recollection for more models and operations conditions


def initialize(values: list, lista: list[str], name: str) -> bool:
    """
    Initialize a condition based on provided parameters.

    The function checks if the provided `name` exists in the `lista`, and if the
    corresponding value from the `values` list indexed by `name` in `lista` equals
    "True". If the condition is satisfied, the function returns True; otherwise,
    returns False.

    Parameters
    ----------
    values : list
        A list of values whose indices correspond to the elements of the `lista`.

    lista : list of str
        A list of string identifiers used to locate the value in `values` that is
        compared to "True".

    name : str
        The specific name or identifier to locate within `lista`.

    Returns
    -------
    bool
        Returns `True` if the condition is satisfied, `False` otherwise.
    """
    if name in lista and values[lista == name] == "True":
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
    :ivar ns_catalyst: Matrix indicating presence of adsorbed species in each catalyst.
    :type ns_catalyst: numpy.ndarray
    :ivar G_formation_rct: Gibbs free energy of formation for reactants.
    :type G_formation_rct: numpy.ndarray
    :ivar G_formation_ads: Gibbs free energy of formation for adsorbed species.
    :type G_formation_ads: numpy.ndarray
    :ivar G_formation_prd: Gibbs free energy of formation for products.
    :type G_formation_prd: numpy.ndarray
    """

    def __init__(self, species_file: str, parameters: object, writer: object) -> None:
        """
        Represents a chemical species processor and organizer for various categories such
        as reactants, products, adsorbed species, and catalysts based on input data
        from a file. This class also processes initial concentrations and thermochemical
        data when applicable.

        Parameters
        ----------
        species_file : str
            Path to the input file containing species data.
        parameters : object
            Object containing simulation or experiment parameters, including flags like
            `thermochemical` and `g_formation` to determine if thermochemical processing
            is required.
        writer : object
            Logger or writer object for real-time messaging and updates during
            processing.
        """
        writer.message(f"Reading species data from file: {species_file}")
        header, raw_data = Collector.raw_data(species_file, writer=writer)
        Collector.column_exists("Species", header, species_file, writer)
        species_list = raw_data[:, header.index("Species")]
        Collector.column_exists("RPACe", header, species_file, writer)
        index = header.index("RPACe")
        self.reactants = species_list[raw_data[:, index] == "R"].tolist()
        self.products = species_list[raw_data[:, index] == "P"].tolist()
        self.adsorbed = species_list[raw_data[:, index] == "A"].tolist()
        self.catalyst = species_list[raw_data[:, index] == "C"].tolist()

        writer.message(
            f"Species lists:\n"
            f"\tReactants: {self.reactants}\n"
            f"\tProducts: {self.products}\n"
            f"\tAdsorbed species: {self.adsorbed}\n"
            f"\tCatalysts: {self.catalyst}"
        )

        Collector.column_exists("c0", header, species_file, writer)
        self.c0_reactants = np.array(
            raw_data[:, header.index("c0")][raw_data[:, index] == "R"].astype(float)
        )
        self.c0_products = np.array(
            raw_data[:, header.index("c0")][raw_data[:, index] == "P"].astype(float)
        )
        self.list = (
                self.reactants + self.products + self.adsorbed + self.catalyst + ["e-"]
        )
        writer.message("Initial concentrations processed.")

        self.ns_catalyst = np.zeros((len(self.catalyst), len(self.adsorbed)))
        if "Sites" in header:
            sites = raw_data[:, header.index("Sites")]
        else:
            sites = np.ones_like(raw_data[:, header.index("c0")], dtype=int)

        for i in range(len(self.catalyst)):
            species_in_catalyst = species_list[
                raw_data[:, header.index("Catalyst")] == self.catalyst[i]
                ]
            nsites = sites[raw_data[:, header.index("Catalyst")] == self.catalyst[i]]
            for specie, ns in zip(species_in_catalyst, nsites):
                if specie in self.adsorbed:
                    self.ns_catalyst[i, self.adsorbed.index(specie)] = float(ns)
        writer.message("Catalyst matrix created.")

        if parameters.thermochemical:
            Collector.column_exists("DG_formation", header, species_file, writer)
            if parameters.g_formation:
                self.g_formation_rct = np.array(
                    raw_data[:, header.index("DG_formation")][
                        raw_data[:, index] == "R"
                        ].astype(float)
                )
                self.g_formation_ads = np.array(
                    raw_data[:, header.index("DG_formation")][
                        raw_data[:, index] == "A"
                        ].astype(float)
                )
                self.g_formation_prd = np.array(
                    raw_data[:, header.index("DG_formation")][
                        raw_data[:, index] == "P"
                        ].astype(float)
                )
            writer.message("Formation energies processed")


class DataReactions:
    """
    Represents a system to process and manage reaction data from a given input file.

    The class extracts and organizes reaction data into matrices for further processing,
    depending on specified parameters and species data. It is designed to handle
    stoichiometric coefficients, catalysts, adsorbates, and other chemical or experimental
    parameters for reaction studies.

    Attributes
    ----------
    list : list
        List of reaction IDs extracted from the input file.
    beta : numpy.ndarray
        Array of Beta values for the reactions, interpreted as floats.
    upsilon : numpy.ndarray
        Reaction matrix including all stoichiometric coefficients and catalysts.
    ne : numpy.ndarray
        Number of electrons transferred for each reaction.
    upsilon_c : numpy.ndarray
        Reaction coefficients matrix excluding catalysts.
    upsilon_a : numpy.ndarray
        Matrix of coefficients for adsorbed species only.
    upsilonx : numpy.ndarray
        Reaction matrix selected for specific model usage (e.g. without catalysts or
        limited to adsorbates), determined by input parameters.
    k_f : numpy.ndarray, optional
        Array of forward kinetic rate constants for experimental data (if applicable).
    k_b : numpy.ndarray, optional
        Array of backward kinetic rate constants for experimental data (if applicable).
    ga : numpy.ndarray, optional
        Free energy of activation values for thermochemical data (if applicable).
    dg_reaction : numpy.ndarray, optional
        Gibbs free energy changes for reactions, if available in data and requested.
    """

    def __init__(
            self, reaction_file: str, parameters: object, species: object, writer: object
    ) -> None:
        writer.message(f"Reading Reactions data from file: {reaction_file}")
        header, raw_data = Collector.raw_data(reaction_file, writer=writer)
        Collector.column_exists("id", header, reaction_file, writer)
        self.list = raw_data[:, header.index("id")].tolist()
        Collector.column_exists("Beta", header, reaction_file, writer)
        self.beta = np.array(raw_data[:, header.index("Beta")].astype(float))

        Collector.column_exists("Reactions", header, reaction_file, writer)
        self.upsilon = np.zeros((len(self.list), len(species.list)))
        for i in range(len(self.list)):
            r = raw_data[:, header.index("Reactions")][i]
            left, right = re.split(r"<->", r)

            species_in_reaction, stoichiometric = self.process_reaction(
                left, species.list
            )
            for specie, coeff in zip(species_in_reaction, stoichiometric):
                self.upsilon[i, species.list.index(specie)] = -coeff

            species_in_reaction, stoichiometric = self.process_reaction(
                right, species.list
            )
            for specie, coeff in zip(species_in_reaction, stoichiometric):
                self.upsilon[i, species.list.index(specie)] = coeff

        self.ne = self.upsilon[:, -1]  # Number of electrons transferred
        self.upsilon = self.upsilon[:, :-1]  # All coefficients, catalysts included
        self.upsilon_c = self.upsilon[
                         :, : -len(species.catalyst)
                         ]  # All coefficients, without catalysts
        self.upsilon_a = self.upsilon_c[
                         :, -len(species.adsorbed):
                         ]  # Adsorbates coefficients
        writer.message("Reaction matrix processed.")

        if parameters.cstr:
            self.upsilonx = self.upsilon_c
        else:
            self.upsilonx = self.upsilon_a

        if parameters.experimental:
            Collector.column_exists("k_f", header, reaction_file, writer)
            Collector.column_exists("k_b", header, reaction_file, writer)
            self.k_f = np.array(raw_data[:, header.index("k_f")].astype(float))
            self.k_b = np.array(raw_data[:, header.index("k_b")].astype(float))
            writer.message("Experimental kinetic rate constants processed.")

        if parameters.thermochemical:
            Collector.column_exists("Ga", header, reaction_file, writer)
            self.ga = np.array(raw_data[:, header.index("Ga")].astype(float))
            if parameters.dg_reaction:
                Collector.column_exists("DG_reaction", header, reaction_file, writer)
                self.dg_reaction = np.array(
                    raw_data[:, header.index("DG_reaction")].astype(float)
                )
            writer.message("Thermochemical reactions parameters processed.")
        # TODO: Add data recollection for more models

    @staticmethod
    def process_reaction(side, species_list):
        """
        Processes a chemical reaction side's string representation and parses the coefficients
        and species in the reaction. Extracted species and their stoichiometric coefficients
        are returned if all species are found in the provided list of known species.

        This method is useful for analyzing and determining which species are involved in
        a chemical reaction and their proportions based on the input reaction string.

        Parameters
        ----------
        side : str
            A string representing one side of a chemical reaction (e.g., '2 H2 + O2').
        species_list : list of str
            A list of valid species names allowed in the chemical reaction.

        Returns
        -------
        tuple of (list of str, list of float)
            A tuple containing two lists:
            - The first list contains the species involved in the reaction as strings.
            - The second list contains their respective stoichiometric coefficients as floats.

        Raises
        ------
        ValueError
            If a species in the reaction string is not present in the species_list or if the
            format of the reaction string cannot be processed.
        """
        species_in_reaction = []
        stoichiometric = []
        str = re.split(r"\s+", side)
        for s in str:
            s = s.strip()
            if s in ["", "+"]:
                continue
            match = re.search(r"(\d+(\.\d+)?)?(.+)", s)
            if match:
                coefficient = float(match.group(1)) if match.group(1) else 1
                specie = match.group(3)
                if specie in species_list:
                    species_in_reaction.append(specie)
                    stoichiometric.append(float(coefficient))
                else:
                    raise ValueError(f"Process reaction error:  {specie} not found")
            else:
                raise ValueError(f"Process reaction error:  {str} not found")
        return species_in_reaction, stoichiometric


class Collector:
    """
    Handles the collection and processing of data files for a project.

    This class serves as the main hub for initializing and managing various data
    components required by the project. It instantiates submodules for handling
    parameters, species, and reactions, while also orchestrating data processing
    from associated files.

    Attributes
    ----------
    directory : str
        The root directory containing the data files and logs.
    parameters : DataParameters
        Manages the parameters extracted from a "parameters.md" file.
    species : DataSpecies
        Manages the species data extracted from a "species.md" file.
    reactions : DataReactions
        Manages the reactions data extracted from a "reactions.md" file.
    """

    def __init__(self, directory, writer: object = None) -> None:
        self.directory = directory
        if writer is None:
            writer = Writer(log_file="melektrodica.log", log_directory=self.directory)
        writer.message("***  Collector  ***")
        self.parameters = DataParameters(
            os.path.join(directory, "parameters.md"), writer
        )
        self.species = DataSpecies(
            os.path.join(directory, "species.md"), self.parameters, writer
        )
        self.reactions = DataReactions(
            os.path.join(directory, "reactions.md"),
            self.parameters,
            self.species,
            writer,
        )
        writer.message("***  Data collection completed successfully.  ***\n")

    @staticmethod
    def raw_data(name_file: str, writer: object):
        """
        Retrieves raw data from a file formatted with a specific delimiter and structure.

        This function reads a file containing data stored in a pipe-separated format. It extracts
        the header from the first row of the file and the data entries starting from the third
        row onwards, while skipping commented or empty lines. The header and data are returned
        separately.

        Parameters
        ----------
        name_file : str
            The file path to be read. It should point to a pipe-separated file containing a
            header and data entries.
        writer : object
            A writer object with a logger attribute capable of logging critical messages.
            The logger is used to log issues related to file accessibility.

        Returns
        -------
        tuple[list[str], numpy.ndarray]
            A tuple containing:
            - A list of strings representing the extracted header.
            - A numpy array containing the raw data entries.

        Raises
        ------
        FileNotFoundError
            Raised if the specified file does not exist at the given file path.
        """
        if not os.path.exists(name_file):
            writer.logger.critical(f"ERROR File not found - {name_file}")
            raise FileNotFoundError(
                f"File was not found at the expected location: {name_file}"
            )
        header = []
        raw_data = []
        with open(name_file, "r") as f:
            lines = f.readlines()
            header = [col.strip() for col in lines[0].split("|")[1:-1]]
            for line in lines[2:]:
                if line.strip() and line[0] != "#":
                    entry_data = [col.strip() for col in line.split("|")[1:-1]]
                    raw_data.append(entry_data)
        raw_data = np.array(raw_data)
        return header, raw_data

    @staticmethod
    def column_exists(
            column_name: str, header: list[str], file_name: str, writer: object
    ):
        """
        Check if a specific column exists within a given header and log an error if it does not.

        This method verifies whether a given column name is present in the provided header.
        If the column is not found, it logs a critical error and raises a ValueError, indicating the absence of
        the required column. The method is static and can be invoked without an instance of the class.

        Parameters
        ----------
        column_name : str
            The name of the column to check for in the header.
        header : list
            A list representing the header of the file, against which the column name is verified.
        file_name : str
            The name of the input file being validated. Used in error logging to specify the source of
            the missing column.
        writer : Any
            A logger or writer object with a `logger.critical()` method for logging critical errors.

        Raises
        ------
        ValueError
            Raised when the required column is not found in the header of the input file.

        """
        if column_name not in header:
            writer.logger.critical(
                f"ERROR: Required column '{column_name}' not found in the header of file '{file_name}'."
            )
            raise ValueError(
                f"Missing required column '{column_name}' in the header of file '{file_name}'."
            )

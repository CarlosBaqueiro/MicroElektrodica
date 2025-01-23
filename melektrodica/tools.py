"""

    μElektrodica© 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Tools

"""

import shutil


class Tool:
    """
    Provides a collection of static utility methods for unit conversion, data visualization,
    text manipulation, displaying introductory banners, LaTeX formatting of chemical species,
    and reading file contents.

    Detailed description of each method is included in their corresponding docstrings.

    Attributes
    ----------
    None
    """

    @staticmethod
    def unit_conversion(variable, value, unit_in, unit_out):
        """
        Converts a given value from one unit to another based on the specified
        physical variable.

        The function currently supports temperature conversion from Celsius to
        Kelvin. Additional conversions for other variables (e.g., length,
        mass, etc.) and their corresponding SI units are to be implemented.

        Parameters
        ----------
        variable : str
            The physical variable for which the conversion is to be performed
            (e.g., "Temperature").
        value : float
            The numerical value to be converted.
        unit_in : str
            The unit of the input value (e.g., "C" for Celsius).
        unit_out : str
            The desired unit for the output value (e.g., "K" for Kelvin).

        Returns
        -------
        float
            The converted value in the desired unit.
        """
        if variable == "Temperature":
            if unit_in == "C" and unit_out == "K":
                value += 273.15
            # TODO: Add conversions for other variables and SI units
        return value

    @staticmethod
    def showme(name, array):
        """
        Displays the name of the input data and details related to the given array such
        as its dimensions and its content.

        This static method is primarily for printing basic array information in a
        readable format, including the name, size (shape), and the entire contents.

        Parameters
        ----------
        name : str
            The name or label to be displayed for the given array. Represents a
            user-defined identification or context of the array's content.
        array : numpy.ndarray
            The input array whose shape (dimensions) and content are to be printed.
            Expected to be a valid NumPy array.

        """
        print("\n", name)
        print("\nSize", array.shape)
        print("\n", array)

    @staticmethod
    def print_center(text):
        """
        Prints the given text centered within the terminal's current width.

        This method calculates the middle position of the terminal width by retrieving
        the number of columns available. Then, it adds a calculated amount of padding
        spaces before the text to align it to the center. Finally, it displays the
        centered text.

        Parameters
        ----------
        text : str
            The text to be printed centered on the terminal.
        """
        columns = shutil.get_terminal_size().columns
        padding = (columns - len(text)) // 2
        text = " " * padding + text
        print(text)

    @staticmethod
    def begin():
        """
        Print an introductory banner for the μElektrodica toolbox.

        This method prints a formatted message containing the name, version, and
        additional information about the μElektrodica electrochemistry toolbox. It is
        designed to display an introductory banner when the toolbox starts.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The banner includes:
        - Toolbox name: μElektrodica
        - Version: 1.0.0
        - Purpose: A Python Electrochemistry Toolbox for Modeling Microkinetic
          Electrocatalytic Reactions
        - Author credits: C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        - License information: CC BY-NC-SA 4.0

        This method relies on the `print_center` method from the `Tool` class to
        display centered text for certain portions of the banner.
        """
        print(f"\n")
        Tool.print_center("\u03BCElektrodica © 2024")
        Tool.print_center("Uxmal 1.0.0\n")
        print(
            f"A Python Electrochemistry Toolbox for Modeling Microkinetic Electrocatalytic Reactions\n"
            f"C. Baqueiro Basto, M. Secanell, L.C. Ordoñez\n"
            f"licensed under CC BY-NC-SA 4.0 \n"
        )

    @staticmethod
    def format_latex_chemical(species):
        """
        Formats a list of chemical species names into LaTeX string representations. The function
        formats individual species to include subscript notation for digits and superscript notation
        for positive charges (e.g., "+"). The resulting formatted names are suitable for inclusion
        in LaTeX documents.

        Parameters
        ----------
        species : list of str
            A list of chemical species represented as strings.

        Returns
        -------
        list of str
            A list of LaTeX-formatted strings where each string corresponds to the
            LaTeX representation of the respective chemical species in the input list.
        """
        chemicals = []
        for chem in species:
            # Añade subíndices: coloca lo que está después de un número como "_" para subíndice
            formatted = ""
            for char in chem:
                if char.isdigit():
                    formatted += f"_{char}"  # Usa subíndice en LaTeX
                elif char == "+":
                    formatted += "^+"
                else:
                    formatted += char
            chemicals.append(
                f"$\\mathrm{{{formatted}}}$"
            )  # Agregar delimitadores de LaTeX
        return chemicals

    @staticmethod
    def show_file(fname):
        """
        Reads the content of a file and returns it as a string.

        This static method opens a file in read mode, reads its entire content, and returns
        it as a string. It is intended for quickly accessing file contents for further
        processing or analysis. The method assumes the file exists and is accessible.

        Parameters
        ----------
        fname : str
            The name of the file to read.

        Returns
        -------
        str
            The contents of the file as a string.
        """
        with open(fname, "r") as f:
            file = f.read()
        return file

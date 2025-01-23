"""

    μElektrodica© 2025
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

        Tools

"""

import shutil


class Tool:
    @staticmethod
    def unit_conversion(variable, value, unit_in, unit_out):
        if variable == "Temperature":
            if unit_in == "C" and unit_out == "K":
                value += 273.15
            # TODO: Add conversions for other variables and SI units
        return value

    @staticmethod
    def showme(name, array):
        print("\n", name)
        print("\nSize", array.shape)
        print("\n", array)

    @staticmethod
    def print_center(text):
        columns = shutil.get_terminal_size().columns
        padding = (columns - len(text)) // 2
        text = " " * padding + text
        print(text)

    @staticmethod
    def begin():
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

    @staticmethod
    def show_file(fname):
        with open(fname, "r") as f:
            file = f.read()
        return file

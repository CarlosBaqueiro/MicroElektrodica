import os
import csv
import logging
from tabulate import tabulate
from colorlog import ColoredFormatter


class Writer:
    _instance = None

    def __new__(cls, log_file="output.log", log_directory="./logs"):
        """
        Singleton pattern implementation to ensure a single instance of the class.
        """
        if cls._instance is None:
            cls._instance = super(Writer, cls).__new__(cls)
            cls._instance._is_initialized = False
        return cls._instance

    def __init__(self, log_file="output.log", log_directory="./logs"):
        """
        Writer class constructor, which initializes the global logger if it is not already configured.
        """
        if not self._is_initialized:
            self.log_file = log_file
            self.log_directory = log_directory
            self.logger = None
            self._setup_logger()
            self._is_initialized = True  # Marks the instance as initialized

            try:
                os.makedirs('MicroElektrodicaResults', exist_ok=True)
            except Exception as e:
                print(f"An error occurred while creating the directory MicroElektrodicaResults: {e}")
                return

    def _setup_logger(self):
        """
        Configures the global logger with support for both console and file logging. Cleans up previous handlers.
        """
        # Set up the main logger
        self.logger = logging.getLogger("GlobalWriterLogger")
        self.logger.setLevel(logging.INFO)

        # Clean up any previous handlers if they exist
        while self.logger.hasHandlers():
            self.logger.removeHandler(self.logger.handlers[0])

        # Ensure the log directory exists
        os.makedirs(self.log_directory, exist_ok=True)
        log_file_path = os.path.join(self.log_directory, self.log_file)

        # Configure the file handler (saves logs to a file)
        file_handler = logging.FileHandler(log_file_path, mode='w')  # Overwrites the file every time
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)

        # Check if running in Jupyter Notebook
        try:
            get_ipython  # This exists only in Jupyter/IPython environments
            is_jupyter = True
        except NameError:
            is_jupyter = False

        # Configure the console handler (logs to the terminal) only if NOT in Jupyter
        if not is_jupyter:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)

            console_formatter = ColoredFormatter(
                "%(log_color)s%(levelname)s: %(message)s",
                log_colors={
                    'DEBUG': 'white',
                    'INFO': 'green',
                    'WARNING': 'yellow',
                    'ERROR': 'red',
                    'CRITICAL': 'bold_red',
                }
            )
            console_handler.setFormatter(console_formatter)
            self.logger.addHandler(console_handler)

    def message(self, message):
        """
        Writes a message to the log file and prints it to the console.
        :param message: Message to write.
        """
        if self.logger:
            self.logger.info(message)

    def to_markdown(self, filname, variable, calculator):
        """
        Write data in a markdown file in table format.
        :param fname: Path to the markdown file (.md) where the data will be written.
        :param variable: The variable type that determines the specific data to write (e.g., "theta" or "j").
        :param results: Object containing the data (e.g., potential, theta, j) to be written into the markdown file.
        """

        fname = os.path.join(calculator.data.directory, 'MicroElektrodicaResults', filname + '.md')
        # Define the columns with "Potential" as the first one
        columns = ["Overpotential [V]"]
        rows = []

        # Define a helper function to format numbers in scientific notation with 3 decimals
        def format_scientific(value):
            return f"{value:.5e}"

        if variable == 'theta' or variable == 'fval':
            # Add column names for adsorbed species
            columns += [specie for specie in calculator.species.adsorbed]
            # Build rows with potential and theta values
            for potential, theta_values in zip(calculator.potential, calculator.results.theta):
                # Convert all values to scientific notation
                formatted_row = [format_scientific(potential)] + [format_scientific(val) for val in theta_values]
                rows.append(formatted_row)
        elif variable == 'j':
            # Add column name for "Current"
            columns += ["Current [A/cm2]"]
            # Build rows with potential and current values
            for potential, current in zip(calculator.potential, calculator.results.j):
                # Convert all values to scientific notation
                formatted_row = [format_scientific(potential), format_scientific(current)]
                rows.append(formatted_row)

        try:
            # Open the markdown file in write mode
            with open(fname, 'w') as md_file:
                # Write the column headers
                md_file.write('| ' + ' | '.join(columns) + ' |\n')
                # Write the separator for the table
                md_file.write('|' + '---|' * len(columns) + '\n')
                # Write the data rows
                for row in rows:
                    md_file.write('| ' + ' | '.join(row) + ' |\n')
            print(f"Successfully written table to {fname}")
        except Exception as e:
            print(f"An error occurred while writing to the file: {e}")

    def to_csv(self, filename, variable, calculator):
        """
        Write data to a CSV file.
        :param filename: Path to the CSV file where the data will be written.
        :param variable: The variable type that determines the specific data to write (e.g., "theta" or "j").
        :param calculator: Object containing the data (e.g., potential, theta, j) to be written into the CSV file.
        """
        fname = os.path.join(calculator.data.directory, 'MicroElektrodicaResults', filename + '.csv')
        # Define the columns with "Potential" as the first one
        columns = ["Overpotential [V]"]
        rows = []

        # Helper function to format numbers in scientific notation with 3 decimals
        def format_scientific(value):
            return f"{value:.5e}"

        # Check which variable type we are handling
        if variable == 'theta' or variable == 'fval':
            # Add column names for adsorbed species
            columns += [specie for specie in calculator.species.adsorbed]
            # Build rows with potential and theta values
            for potential, theta_values in zip(calculator.potential, calculator.results.theta):
                # Convert all values to scientific notation
                formatted_row = [format_scientific(potential)] + [format_scientific(val) for val in theta_values]
                rows.append(formatted_row)
        elif variable == 'j':
            # Add column name for "Current"
            columns += ["Current [A/cm2]"]
            # Build rows with potential and current values
            for potential, current in zip(calculator.potential, calculator.results.j):
                # Convert all values to scientific notation
                formatted_row = [format_scientific(potential), format_scientific(current)]
                rows.append(formatted_row)

        try:
            # Open the CSV file in write mode
            with open(fname, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file)
                # Write the column headers
                writer.writerow(columns)
                # Write the data rows
                writer.writerows(rows)
            print(f"Successfully written data to {fname}")
        except Exception as e:
            print(f"An error occurred while writing to the file: {e}")

    def display_table(self, variable, fitter):
        """
        Displays a table with names and their associated energies, based on the selected variable type.

        :param variable: The type of data to display ('Reactions' or 'Species').
        :param fitter: An object containing the necessary data (reactions, species, and energy values).
        :raises ValueError: If the 'variable' parameter is invalid.
        """
        if variable == 'Reactions':
            items = fitter.data.reactions.list
            energies = fitter.ga_fit
            col = 'Free energies of activation reactions [eV]'
        elif variable == 'Species':
            items = fitter.data.species.adsorbed
            energies = fitter.gf_fit
            col = 'Formation energy [eV]'
        else:
            raise ValueError("The 'variable' parameter must be 'Reactions' or 'Species'.")

        # Define table columns and rows
        columns = [variable, col]
        rows = [[item, energy] for item, energy in zip(items, energies)]

        return tabulate(rows, headers=columns, tablefmt="grid")
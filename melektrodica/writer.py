import logging
import os
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

        # Configure the console handler (logs to the terminal)
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

    def markdown(self, filname, variable, calculator):
        """
        Write data in a markdown file in table format.
        :param fname: Path to the markdown file (.md) where the data will be written.
        :param variable: The variable type that determines the specific data to write (e.g., "theta" or "j").
        :param results: Object containing the data (e.g., potential, theta, j) to be written into the markdown file.
        """

        fname = os.path.join(calculator.data.directory, filname + '.md')
        # Define the columns with "Potential" as the first one
        columns = ["Potential (V)"]
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
            columns += ["Current (A/cm2)"]
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
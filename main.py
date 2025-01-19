"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

"""

import os
import argparse
from melektrodica import *

Tool.begin()


def default_calculations(directory):
    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
    else:
        writer = Writer(log_file="melektrodica.log", log_directory=directory)
        print(f"Directory: {directory}")
        data = Collector(directory)


def adjust_curves(directory, file):
    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
        pass


if __name__ == "__main__":
        # Create the parser for command-line arguments
    parser = argparse.ArgumentParser(description="Script to perform calculations and curve fitting.")

    # Positional argument: directory
    parser.add_argument("directory", type=str, help="Main directory to process the data.")

    # Optional argument: data file for curve fitting
    parser.add_argument("-f", "--file", type=str, help="Data file for curve fitting.")

    # Parse arguments
    args = parser.parse_args()

    # Logic for handling the arguments
    if args.file:
        # Case 2: Curve fitting
        print(f"Curve fitting enabled. Processing the data file: {args.file}")
        # Add your code for curve fitting using args.file
        adjust_curves(args.directory, args.file)
    else:
        # Case 1: Default calculations
        print(f"Performing default calculations in the directory: {args.directory}")
        # Add your code for default calculations
        default_calculations(args.directory)

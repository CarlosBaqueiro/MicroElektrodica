"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

"""


import os
import argparse
from Source import *
begin()

def run(directory):
    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
    else:
        print(f"Directory: {directory}")
        data = Collector(directory)
        Grapher(data, Calculator(data))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str, help='Directory with the input files')
    args = parser.parse_args()
    run(args.directory)

"""

    μElektrodica (Uxmal, version 1.0.0)
        A Python Tool for Modeling Microkinetic Electrocatalytic Reactions
        Copyright (C) 2024 C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        All rights reserved.

"""


import os
from Source import *
begin()

pakage = 'Wang2007Hydrogen'

run = True
validation = False

if __name__ == '__main__':

    if run:
        directory = os.path.join('Examples', pakage)
        if not os.path.exists(directory):
            print(f"The directory {directory} does not exist.")
        else:
            print(f"Directory: {directory}")
            data = Collector(directory)
            Grapher(data, Calculator(data))
    if validation:
        Validation()



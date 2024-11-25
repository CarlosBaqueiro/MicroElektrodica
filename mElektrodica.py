"""

    μElektrodica © 2024
        by C. Baqueiro Basto, M. Secanell, L.C. Ordoñez
        is licensed under CC BY-NC-SA 4.0

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



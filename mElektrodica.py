"""
μElektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import os
from Source import *
begin()

pakage = 'Wang2007Hydrogen'
#pakage = 'SanchezMonreal2017Ethanol'
#pakage =  'Moore2013Oxygen'
if __name__ == '__main__':
    Validation()
    run = False
    test = False
    if run:
        directory = os.path.join('Examples', pakage)
        if not os.path.exists(directory):
            print(f"The directory {directory} does not exist.")
        else:
            print(f"Directory: {directory}")
            data = Collector(directory)
            results = Calculator(data)
            Grapher(data, results)




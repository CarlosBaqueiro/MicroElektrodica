import shutil

def unit_conversion(variable, value, unit_in, unit_out):
    if variable == 'Temperature':
        if unit_in == 'C' and unit_out == 'K':
            value += 273.15
        #TODO: Add conversions for other variables and SI units
    return value

def showme (name, array):
    print("\n", name)
    print("\nSize", array.shape)
    print("\n", array)

def print_center(text):
    columns = shutil.get_terminal_size().columns
    padding = (columns - len(text)) // 2
    text = ' ' * padding + text
    print(text)

def begin ( ):
    print(f'\n')
    print_center('\u03BCElektrodica')
    print_center('Uxmal 1.0.0\n')
    print(f'A Python Tool for Modeling Microkinetic Electrocatalytic Reactions\n'
          f'Copyright (C) 2024 C. Baqueiro Basto, M. Secanell, L.C. Ordoñez\n'
          f'All rights reserved.\n')
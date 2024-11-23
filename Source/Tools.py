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
    print(f'\n')
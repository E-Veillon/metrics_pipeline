#!/usr/bin/python

'''Counts the number of structures in a given CIF file.'''

from argparse import ArgumentParser
from pathlib import Path

def has_str(string: str, pattern: str) -> bool:
    '''
    Returns True if pattern is in string, otherwise returns False

    Parameters:
        string (str): string to search in.
        pattern (str): which string to search.
    
    Returns:
        Whether the searched string was found or not.
    '''
    find_index = string.find(pattern)
    return not (find_index == -1)

def main():
    parser = ArgumentParser(
    prog='cif_counter.py', 
    description='Counts the number of structures in a given CIF file.'
    )

    parser.add_argument(
    'filename', 
    type=str, 
    help='file to count structures in.'
    )

    args = parser.parse_args()
    file = args.filename

    assert Path(file).is_file() and file.endswith('.cif'), \
    'Provided filename must be a valid file in CIF format.'

    with open(file,'r') as data:
        nbr_structs = len(list(filter(lambda x: has_str(x, 'data_'), list(data))))

    print(f"{nbr_structs} structures found in file '{file}'.")


if __name__ == '__main__':
    main()
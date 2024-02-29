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


    # ARGUMENTS PARSING BLOCK

    parser = ArgumentParser(
    prog='cif_counter.py', 
    description='Counts the number of structures in a given CIF file.'
    )

    parser.add_argument(
    'filename', 
    type=str, 
    help='file to count structures in.'
    )
    parser.add_argument(
    '--cut',
    type=int,
    help='create a truncated copy of the tested file containing only given number of structures, going from the first one.',
    metavar='int'
    )

    args = parser.parse_args()
    filename = args.filename
    cut_nbr  = args.cut if isinstance(args.cut, int) and (args.cut >= 1) else None

    assert Path(filename).is_file() and filename.endswith('.cif'), \
    'Provided filename must be a valid file in CIF format.'


    # MAIN BLOCK

    with open(filename,'r') as data:
        lines = data.read().splitlines()

    data_breakpoints = list(filter(lambda x: has_str(x, 'data_'), lines))
    nbr_structs = len(data_breakpoints)

    if cut_nbr is not None and cut_nbr >= nbr_structs:
        print(
            f"Provided 'cut' arg ({cut_nbr}) is larger or equal to the total number of structures in {filename}.\n"
            "Therefore, as it will not change anything, the cut option will be deactivated."
        )

    elif cut_nbr is not None: 
        wrong_cut = False
        idx_changed_lines = []
        cut_line = data_breakpoints[cut_nbr]

        if data_breakpoints.index(cut_line) < cut_nbr: # Une ligne identique apparaît avant celle voulue
            nbr_iter = 0
            while True:
                idx_same_line = data_breakpoints.index(cut_line)
                if idx_same_line == cut_nbr:
                    break
                elif idx_same_line < cut_nbr:
                    idx_changed_lines.append(idx_same_line)
                    data_breakpoints[idx_same_line] = f'data_xx{nbr_iter}xx'
                    lines[lines.index(cut_line)] = data_breakpoints[idx_same_line]
                    nbr_iter += 1
                    continue
                else:
                    print(
                        'Unexpected behaviour happened while trying to cut the file.\n'
                        'The cutting option will stop here without doing anything.'
                    )
                    wrong_cut = True
                    break

        if not wrong_cut:
            cut_lines = lines[:lines.index(cut_line)]

            for idx in idx_changed_lines:
                line_to_restore = data_breakpoints[idx]
                idx_to_restore = lines.index(line_to_restore)
                lines[idx_to_restore] = cut_line

            cut_text = '\n'.join(cut_lines)
            cut_file = filename.replace('.cif', f'_{cut_nbr}.cif')

            with open(cut_file, 'wt') as outfile:
                outfile.write(cut_text)

            print(f"The file '{cut_file}' containing the first {cut_nbr} structure(s) from '{filename}' was successfully created.")

    print(f"{nbr_structs} structures found in file '{filename}'.")


if __name__ == '__main__':
    main()

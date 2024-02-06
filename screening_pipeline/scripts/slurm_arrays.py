#!/usr/bin/python
'''
Small IO script to test Slurm job arrays.
'''

import os
from argparse import ArgumentParser, RawTextHelpFormatter
from datetime import datetime
from pathlib import Path


def main():
    start = datetime.now()

    parser = ArgumentParser(prog='slurm_arrays.py', formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        'input', 
        type=str, 
        help='Input file.'
    )
    parser.add_argument(
        'index', 
        type=int, 
        help='Index number of the line to write. The output file will have the form "_[index_nbr].txt".\n \
            This argument is meant to be used with slurm job array tasks indexes.', 
        metavar='job_array_task_id' 
    )
    parser.add_argument(
        '-o', 
        '--outdir', 
        type=str, 
        default='[input]_out', 
        help='Optional output directory. Defaults to same name as input file with "_out" added at the end.'
    )

    args = parser.parse_args()

    assert Path(args.input).is_file(), \
    f'{args.input}: No such file found.'

    assert args.input.endswith('.txt'), \
    f'Input file must be a .txt file (got {args.input})'

    assert args.index >= 0, \
    f'Line index must be positive or zero.'

    if args.outdir == '[input]_out':
        args.outdir = args.input.replace('.txt', '_out')
    
    Path(args.outdir).mkdir(exist_ok=True)

    with open(args.input, 'rt') as read_file:
        lines = read_file.readlines()

    assert args.index < len(lines), \
    f'Provided index "{args.index}" is out of the range of the text ({len(lines)-1}).'

    with open(os.path.join(args.outdir, f'_{args.index}.txt'), 'wt') as out_file:
        out_file.write(lines[args.index])
    
    stop = datetime.now()
    print(f'Elapsed time: {stop-start}')


if __name__=='__main__':
    main()

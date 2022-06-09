#!/usr/bin/env python3
"""
Finds cliquemates or grammates.
"""

__author__ = 'Yosuke Mizutani'
__version__ = '0.0.1'
__license__ = 'Apache License, Version 2.0'

# imports standard libraries
import sys
import argparse
import mates_finder_lib as lib


def get_parser():
    """Argument parser."""

    parser = argparse.ArgumentParser(description='Finds cliquemates or grammates.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('--type', choices=['clique', 'line', 'gram'], default='gram', help='type of the mates to find')
    parser.add_argument('--n', type=int, default=2, help='number of vertices')
    parser.add_argument('--m', type=int, default=2, help='number of hyperedges')

    return parser


def main(args):
    """Entry point of the program. """

    # main logic
    if args.type == 'gram':
        print(f'[Grammates: n={args.n}, m={args.m}]')
        lib.print_nontrivial_grammates(args.n, args.m)
        print('Done')
    elif args.type == 'clique':
        print(f'[Cliquemates: n={args.n}, max_m={args.m}]')
        ret = lib.get_largest_cliquemates(args.n, args.m)
        print(f'size: {len(ret)}')
        print(ret)
    elif args.type == 'line':
        print(f'[Linemates: max_n={args.n}, m={args.m}]')
        ret = lib.get_largest_linemates(args.n, args.m)
        print(f'size: {len(ret)}')
        print(ret)


if __name__ == '__main__':
    main(get_parser().parse_args())

#!/usr/bin/env python

"""
Split a VCF file into N equally sized batches
"""

import argparse
import os
import sys
import gzip
from enum import Enum

from somatypus_utility_lib import SomatypusError

class Strategy(Enum):
    Rows = 'Rows'
    Files = 'Files'

def parse_args():
    parser = argparse.ArgumentParser(description=('Split a VCF file into N'
                                                  ' equally sized batches'))
    parser.add_argument('--vcf', type=str, required=True
                        , help='VCF file to split')
    parser.add_argument('-n', '--number', type=int, default=10,
                        help='Number of output files')
    parser.add_argument('-o', '--outprefix', type=str, default='out',
                        help='Output file prefix')
    parser.add_argument('-r', '--rows-per-file', type=int, default=2000000,
                        help='Number of rows per output file')
    parser.add_argument('-s', '--strategy', type=Strategy,
                        choices=list(Strategy), default=Strategy.Rows,
                        help=('Strategy to use for splitting the file:' 
                              'Rows (default): split into files with '
                              '`--rows-per-file` rows;'
                              ' Files: split into `--number` files'))
    return parser.parse_args()

def is_gzip(filename):
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def validate_args(args):
    if args.number < 1:
        raise SomatypusError('Number of output files must be greater than 0')

    if not os.path.isfile(args.vcf):
        raise SomatypusError('VCF file does not exist: {}'.format(args.vcf))

def process_file_strategy_n(f, number, outprefix):
    print("Strategy = Files")
    print("Number of files = {}".format(number))
    header_lines = []
    filenames = ['{}-{}.vcf'.format(outprefix, i) for i in range(number)]
    outfiles = [open(outfile, 'w') for outfile in filenames]
    linenum = 0
    try:
        written_header = False
        for linenum, line in enumerate(f):
            if line.startswith('#'):
                header_lines.append(line)
            else:
                if not written_header:
                    header = ''.join(header_lines)
                    for outfile in outfiles:
                        outfile.write(header)
                    written_header = True
                outfiles[linenum % number].write(line)
                if (linenum + 1) % 100000 == 0:
                    sys.stderr.write('Processed {} lines\r'.format(linenum + 1))
                    sys.stderr.flush()
        sys.stderr.write('Processed {} lines\n'.format(linenum + 1))
        return filenames
    except Exception as e:
        raise SomatypusError('Error processing line {}: {}'.format(linenum + 1, e))
    finally:
        for outfile in outfiles:
            outfile.close()

def process_file_strategy_r(f, rows_per_file, outprefix):
    print("Strategy = Rows")
    def new_filename(outprefix, file_index):
        return '{}-{}.vcf'.format(outprefix, file_index) 

    header_lines = []
    file_index = 0
    outfile = None
    rows_written_to_current_file = 0
    rows_written_total = 0
    filenames = []
    try:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                if filenames == [] or rows_written_to_current_file >= rows_per_file:
                    if outfile is not None:
                        outfile.close()
                    filename = new_filename(outprefix, file_index)
                    outfile = open(filename, 'w')
                    outfile.write(''.join(header_lines))
                    filenames.append(filename)
                    file_index += 1
                    rows_written_to_current_file = 0
                outfile.write(line)
                rows_written_to_current_file += 1
                rows_written_total += 1
                if rows_written_total % 100000 == 0:
                    sys.stderr.write('Processed {} lines\r'.format(rows_written_total))
                    sys.stderr.flush()
        outfile.close()
        sys.stderr.write('Processed {} lines\n'.format(rows_written_total))
        return filenames
    except Exception as e:
        raise SomatypusError('Error processing line {}: {}'.format(rows_written_total + 1, e))
    finally:
        if not outfile.closed:
            outfile.close()

def main():
    args = parse_args()
    try:
        validate_args(args)
        filenames = []
        if is_gzip(args.vcf):
            with gzip.open(args.vcf, 'rt') as f:
                if args.strategy == Strategy.Files:
                    filenames = process_file_strategy_n(f, args.number, args.outprefix)
                else:
                    filenames = process_file_strategy_r(f, args.rows_per_file, args.outprefix)
        else:
            with open(args.vcf, 'r') as f:
                if args.strategy == Strategy.Files:
                    filenames = process_file_strategy_n(f, args.number, args.outprefix)
                else:
                    filenames = process_file_strategy_r(f, args.rows_per_file, args.outprefix)
        if is_gzip(args.vcf):
            for f in filenames:
                sys.stderr.write('Compressing {}...\n'.format(f))
                os.system('bgzip {}'.format(f))
                os.system('tabix --csi {}.gz'.format(f))
        sys.stderr.write('Done\n')
    except SomatypusError as e:
        sys.stderr.write('ERROR: {}\n'.format(e))
        sys.exit(1)
    except Exception as e:
        sys.stderr.write('UNEXPECTED FATAL ERROR: {}\n'.format(e))
        sys.exit(2)

if __name__ == '__main__':
    main()

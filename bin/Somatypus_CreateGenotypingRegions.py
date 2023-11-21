#!/usr/bin/env python

import gzip
import os
import sys
from collections import OrderedDict

class ProgramArgException(Exception):
    pass

class ProgramLogicException(Exception):
    pass

def is_gzip(filename):
    """
    Check the first two bytes of the file to see if
    it is gzipped
    """
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

class FileReader:
    def __init__(self, filename):
        if is_gzip(filename):
            self._file = gzip.open(filename, 'rt')
        else:
            self._file = open(filename, 'r')

    def __enter__(self):
        return self._file

    def __exit__(self, type, value, traceback):
        if not self._file.closed:
            self._file.close()

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Create genotyping regions')
    parser.add_argument('filename', help='VCF/VCF.gz file')
    parser.add_argument('-o', '--output', help='Output file name', default='regions.txt')
    parser.add_argument('-w', '--window', type=int, default=100000, help='Window size')
    return parser.parse_args()

def check_args(args):
    if args.window < 1:
        raise ProgramArgException('Window size must be greater than 0')
    if not os.path.isfile(args.filename):
        raise ProgramArgException('File {} does not exist'.format(args.filename))

def print_user_info(args):
    sys.stdout.write('Creating genotyping regions for {}\n'.format(args.filename))
    sys.stdout.write('Using window size = {}\n'.format(args.window))
    sys.stdout.write('Output will go to {}\n'.format(args.output))

def run(args):
    regions = OrderedDict()
    with FileReader(args.filename) as fl:
        for line in fl:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            start = pos - (pos % args.window) + 1
            end = start + args.window - 1
            if end - start + 1 != args.window:
                raise ProgramLogicException(
                        'Window size is not right: {}:{} start={} end={} window={}'.format(
                            chrom, pos, start, end, args.window))
            if pos < start or pos > end:
                raise ProgramLogicException(
                        'Window doesn\'t enclose position: {}:{} start={} end={} window={}'.format(
                            chrom, pos, start, end, args.window))

            regions[(chrom, start, end)] = pos
    return regions

def write_output(args, regions):
    with open(args.output, 'w') as fl:
        for (chrom, start, end) in regions.keys():
            fl.write('{}:{}-{}\n'.format(chrom, start, end))

def main():
    args = parse_args()
    try:
        check_args(args)
        print_user_info(args)
        regions = run(args)
        write_output(args, regions)
        sys.stdout.write('Done\n')

    except ProgramArgException as e:
        sys.stderr.write('CommandLine Error: {}\n'.format(e))
        sys.exit(1)

    except ProgramLogicException as e:
        sys.stderr.write('Program Error: {}\n'.format(e))
        sys.exit(2)

    except Exception as e:
        sys.stderr.write('Unexpected Error: {}\n'.format(e))
        sys.exit(3)

if __name__ == '__main__':
    main()

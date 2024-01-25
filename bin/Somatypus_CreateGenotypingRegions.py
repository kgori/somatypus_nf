#!/usr/bin/env python

"""
Creates genotyping regions for Platypus from a VCF file

The output file will contain a list of regions in the format:
   chr:start-end

Also,
1: Fits the region to the size of the variant (i.e. indel size, or 2bp for doublet SNPs)
2: Adds a window of '-w' bases around each variant
3: Merges overlapping windows
4: Sorts the windows by chromosome and position
"""
import gzip
import os
import sys
from somatypus_utility_lib import Interval, get_variant_from_line

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
    parser.add_argument('-w', '--window', type=int, default=1000, help='Window size')
    return parser.parse_args()

def check_args(args):
    if args.window < 0:
        raise ProgramArgException('Window size must be 0 or greater')
    if not os.path.isfile(args.filename):
        raise ProgramArgException('File {} does not exist'.format(args.filename))

def print_user_info(args):
    sys.stdout.write('Creating genotyping regions for {}\n'.format(args.filename))
    sys.stdout.write('Using window size = {}\n'.format(args.window))
    sys.stdout.write('Output will go to {}\n'.format(args.output))

def abs(x):
    return x if x >= 0 else -x

def run(args):
    chrs = []
    regions = {}

    with FileReader(args.filename) as fl:
        for line in fl:
            if line.startswith('#'):
                continue
            variant = get_variant_from_line(line)

            if variant.is_multiallelic():
                raise ProgramLogicException('Multiallelic variants are not supported')

            if variant.chrom not in chrs:
                chrs.append(variant.chrom)
                regions[variant.chrom] = []

            variant_size = abs(len(variant.ref) - len(variant.alt)) + 1
            ivl_start = variant.pos
            ivl_end = variant.pos + variant_size - 1

            # Adjust the interval with the window size
            ivl_start = max(1, ivl_start - args.window)
            ivl_end = max(ivl_start, ivl_end + args.window)
            ivl = Interval(ivl_start, ivl_end)

            if not variant.pos in ivl:
                raise ProgramLogicException('Created interval does not contain the variant position: pos={}, ivl={}'.format(variant.pos, ivl))

            if len(regions[variant.chrom]) == 0:
                regions[variant.chrom].append(ivl)
            else:
                prev = regions[variant.chrom][-1]
                if prev.intersects(ivl):
                    prev.absorb(ivl)
                else:
                    regions[variant.chrom].append(ivl)
            # The current pos should be within the most recently added interval
            if not variant.pos in regions[variant.chrom][-1]:
                raise ProgramLogicException('Variant position {} not in {}'.format(variant.pos, regions[variant.chrom][-1]))

    result = []
    for chr in chrs:
        # For each chromosome the intervals should be sorted and disjoint
        if len(regions[chr]) > 1:
            for i in range(len(regions[chr])-1):
                if not regions[chr][i] < regions[chr][i+1]:
                    raise ProgramLogicException('Intervals are unsorted')
                if regions[chr][i].intersects(regions[chr][i+1]):
                    raise ProgramLogicException('Intervals are not disjoint')

        for ivl in regions[chr]:
            result.append((chr, ivl.lower, ivl.upper))

    return result

def write_output(args, regions):
    with open(args.output, 'w') as fl:
        for (chrom, start, end) in regions:
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

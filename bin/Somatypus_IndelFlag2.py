#!/usr/bin/env python
import argparse
import os
import sys
from collections import defaultdict
from enum import Enum

from somatypus_utility_lib import \
    SomatypusError, Interval, Variant

class Sorting(Enum):
    Positional = 'Positional'
    Lexicographical = 'Lexicographical'

    def __str__(self):
        return self.value

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract SNVs near indels from each sample.')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('-w', '--window', type=int, default=5,
                        help='Window size (default: 5)')
    parser.add_argument('-s', '--sorting', type=Sorting, choices=list(Sorting),
                        default=Sorting.Lexicographical,
                        help='Sorting of the output file (default: Lexicographical)')
    return parser.parse_args()

def validate_args(args):
    if not os.path.exists(args.vcf_file):
        raise SomatypusError('VCF file not found: {}\n'.format(args.vcf_file))

def main():
    args = parse_args()

    try:
        validate_args(args)
        intervals = get_indel_intervals(args.vcf_file, args.window)
        variants = get_snps_near_indels(args.vcf_file, intervals)
        if args.sorting == Sorting.Lexicographical:
            variants = sorted(variants, key=lambda v: str(v))
        with open(args.output_file, 'w') as out:
            for variant in variants:
                out.write('{}\n'.format(variant))
    except SomatypusError as e:
        sys.stderr.write('ERROR: {}\n'.format(e))
        sys.exit(1)
    except Exception as e:
        sys.stderr.write('UNEXPECTED FATAL ERROR: {}\n'.format(e))
        sys.exit(2)

def get_variant_from_line(line):
    col = line.strip().split('\t')
    chrom = col[0]
    pos = col[1]
    ref = col[3]
    alt = col[4]
    return Variant(chrom, pos, ref, alt)

def get_indel_intervals(vcf_file, window=5):
    # Extract SNVs near indels from each sample
    sys.stdout.write('[indel_intervals]: Searching file {}\n'.format(vcf_file))
    
    intervals = defaultdict(list) # {chrom: [Interval, Interval, ...]}
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                var = get_variant_from_line(line)
                if var.is_multiallelic():
                    raise SomatypusError(
                            'Multiallelic variant found at {0}:{1}. Use splitMAandMNPs.py first.'
                            .format(var.chrom, var.pos))
                # If indel, mark a window of WINDOW bp at both flanks 
                # of the indel position/footprint
                if var.is_indel():
                    # If deletion: footprint is the length of the REF,
                    # only downstream
                    ftprint = len(var.ref) - 1 if var.is_del() else 0
                    interval = Interval(var.pos - window,
                                        var.pos + ftprint + window)
                    ivl_list = intervals[var.chrom]
                    ivl_last = ivl_list[-1] if ivl_list else None
                    if ivl_last and ivl_last.intersects(interval):
                        ivl_last.absorb(interval)
                    else:
                        ivl_list.append(interval)
    # Ensure intervals are sorted even if vcf file is unsorted
    for k, v in intervals.items():
        intervals[k] = sorted(v)

    return intervals

def get_snps_near_indels(vcf_file, intervals):
    """
    Extract SNVs near indels from each sample. Relies on the VCF being
    sorted by position within each chromosome. If a variant is found out
    of sorted order, the function falls back to a slower algorithm.
    """
    sys.stdout.write('[snps_near_indels]: Searching file {}\n'.format(vcf_file))
    chrom = None
    ivl_list = []
    ivl_iter = None
    ivl = None
    vars = set()
    prev_pos = 0 # check vcf_file is strictly sorted
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                var = get_variant_from_line(line)
                if var.chrom != chrom:
                    chrom = var.chrom
                    prev_pos = 0
                    ivl_list = intervals[chrom]
                    ivl_iter = iter(ivl_list) if ivl_list else None
                    ivl = next(ivl_iter) if ivl_iter else None
                if len(ivl_list) > 0:
                    if var.is_snv():
                        if var.pos < prev_pos:
                            sys.stderr.write('WARNING: VCF file is not strictly sorted. Falling back to slower algorithm for unsorted input.\n')
                            return get_snps_near_indels_unsorted(vcf_file,
                                                                 intervals)
                        prev_pos = var.pos
                        while var.pos > ivl.upper:
                            try:
                                ivl = next(ivl_iter)
                            except StopIteration:
                                break
                        if var.pos in ivl:
                            vars.add(var)
    return sorted(vars)

                    
def get_snps_near_indels_unsorted(vcf_file, intervals):
    """
    Alternative implementation for unsorted VCF files.
    (around 2.5x slower)
    """
    sys.stdout.write('[snps_near_indels_unsorted]: Searching file {}\n'.format(vcf_file))
    from bisect import bisect
    chrom = None
    ivl_list = []
    vars = set()
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                var = get_variant_from_line(line)
                if var.chrom != chrom:
                    chrom = var.chrom
                    ivl_list = intervals[chrom]
                if len(ivl_list) > 0:
                    if var.is_snv():
                        search_result = bisect(ivl_list,
                                           Interval(var.pos, sys.maxsize))
                        ivl_index = max(0, search_result - 1)
                        ivl = ivl_list[ivl_index]
                        if var.pos in ivl:
                            vars.add(var)
    return sorted(vars)


if __name__ == '__main__':
    main()

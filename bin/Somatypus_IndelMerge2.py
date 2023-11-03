#!/usr/bin/env python
import argparse
import os
import sys
from collections import defaultdict


from somatypus_utility_lib import \
    SomatypusError, Interval, Variant, \
    get_variant_from_line


def minimal_vcf_string(variant):
    return '{}\t{}\t.\t{}\t{}\t.\t.\t.'.format(variant.chrom, variant.pos, variant.ref, variant.alt)


def minimal_vcf_header():
    return '##fileformat=VCFv4.1\n##source=Somatypus\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge indels from multiple VCFs')
    subparsers = parser.add_subparsers(dest='command')
    parser_extract = subparsers.add_parser('extract', help='Extract indels from single VCF file')
    parser_extract.add_argument('vcf_file', help='VCF file')
    parser_extract.add_argument('output_file', help='Output file')
    parser_merge = subparsers.add_parser('merge', help='Merge indels from results of extract')
    parser_merge.add_argument('merge_files', nargs='+', help='Files from extract to merge')
    parser_merge.add_argument('output_file', help='Output file')
    return parser.parse_args()

def validate_args(args):
    if args.command == 'extract':
        if not os.path.isfile(args.vcf_file):
            raise SomatypusError('VCF file not found: {}'.format(args.vcf_file))
    elif args.command == 'merge':
        for extract_file in args.merge_files:
            if not os.path.isfile(extract_file):
                raise SomatypusError('Extract file not found: {}'.format(extract_file))
        if os.path.isfile(args.output_file):
            raise SomatypusError('Output file already exists: {}'.format(args.output_file))

def extract_indels(vcf_file):
    indels = set()
    flags = ['badReads', 'MQ', 'strandBias', 'SC', 'QD']
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            variant = get_variant_from_line(line)
            filter = line.strip().split('\t')[6]
            if variant.is_indel() and not variant.is_multiallelic():
                if not any(flag in filter for flag in flags):
                    indels.add(variant)
    return indels

def main():
    args = parse_args()

    try:
        validate_args(args)

        if args.command == 'extract':
            indels = extract_indels(args.vcf_file)
            with open(args.output_file, 'w') as output:
                for indel in sorted(indels):
                    output.write(minimal_vcf_string(indel) + '\n')

        elif args.command == 'merge':
            d = defaultdict(set)
            for filename in args.merge_files:
                print('Processing {}'.format(filename))
                with open(filename) as f:
                    for line in f:
                        variant = get_variant_from_line(line)
                        d[(variant.chrom, variant.pos)].add(variant)
            print('Writing output to {}'.format(args.output_file))
            with open(args.output_file, 'w') as output:
                output.write(minimal_vcf_header() + '\n')
                for _, variants in sorted(d.items()):
                    if len(variants) == 1:
                        variant = variants.pop()
                        output.write(minimal_vcf_string(variant) + '\n')

    except SomatypusError as e:
        sys.stderr.write('ERROR: {}\n'.format(e))
        sys.exit(1)
    except Exception as e:
        sys.stderr.write('UNEXPECTED FATAL ERROR: {}\n'.format(e))
        sys.exit(2)


if __name__ == '__main__':
    main()


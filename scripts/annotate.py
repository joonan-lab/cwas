#!/usr/bin/env python
"""
Script for genomic and functional annotations using Variant Effect Predictor (VEP) and user-added custom BED files
"""

import argparse
import multiprocessing as mp
import os
import yaml

import pysam

from utils import get_curr_time, execute_cmd, bgzip_tabix
from collections import deque


def main():
    print(__doc__)

    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()
    print_args(args)
    check_args_validity(args)
    print()

    # Path settings
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    vep_custom_conf_path = os.path.join(project_dir, 'conf', 'vep_custom.yaml')
    annot_bed_path = os.path.join(project_dir, 'data', 'annotate', 'merged_annotation.bed.gz')
    tmp_vcf_path = args.out_vcf_path.replace('.vcf', '.tmp.vcf')  # Temporary file for a result of VEP
    tmp_vcf_gz_path = tmp_vcf_path + '.gz'
    tmp_vcf_gz_idx_path = tmp_vcf_gz_path + '.tbi'

    # Annotate by Ensembl Variant Effect Predictor (VEP)
    vep_custom_path_dict = {}

    with open(vep_custom_conf_path) as vep_custom_conf_file:
        vep_custom_conf = yaml.safe_load(vep_custom_conf_file)

    for file_key in vep_custom_conf:
        vep_custom_path_dict[file_key] = os.path.join(project_dir, vep_custom_conf[file_key])

    print(f'[{get_curr_time()}, Progress] Run Variant Effect Predictor (VEP)')
    if os.path.isfile(tmp_vcf_gz_path):
        print(f'[{get_curr_time()}, Progress] The temporary VEP result already exists so skip this VEP step')
    else:
        cmd = make_vep_cmd(args.vep_script, args.in_vcf_path, tmp_vcf_path, vep_custom_path_dict)
        execute_cmd(cmd, True)
        bgzip_tabix(tmp_vcf_path)

    # Annotate by the early prepared BED file with merged annotation information
    print(f'[{get_curr_time()}, Progress] Annotate by user-added BED files')
    annotate_by_bed(tmp_vcf_gz_path, args.out_vcf_path, annot_bed_path)

    # Remove the temporary files
    os.remove(tmp_vcf_gz_path)
    os.remove(tmp_vcf_gz_idx_path)

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True, type=str, help='Input VCF file')
    parser.add_argument('-o', '--outfile', dest='out_vcf_path', required=False, type=str,
                        help='Path of the VCF output', default='annotate_output.vcf')
    parser.add_argument('-s', '--split', dest='split_vcf', required=False, type=int, choices={0, 1},
                        help='Split the input VCF by chromosome and run VEP for each split VCF (Default: 0)',
                        default=0)
    parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                        help='Number of processes for this script (only necessary for split VCF files) (Default: 1)',
                        default=1)
    parser.add_argument('--vep', dest='vep_script', required=False, type=str,
                        help='Path of a Perl script to execute VEP (Default: vep (binary))', default='vep')

    return parser


def print_args(args: argparse.Namespace):
    print(f'[Setting] The input VCF file: {args.in_vcf_path}')
    print(f'[Setting] The output path: {args.out_vcf_path}')
    print(f'[Setting] VEP script: {args.vep_script}')

    if args.split_vcf:
        print(f'[Setting] Split the input VCF file by chromosome and run vep for each split VCF')
        print(f'[Setting] No. processes for this script: {args.num_proc:,d}')


def check_args_validity(args: argparse.Namespace):
    assert os.path.isfile(args.in_vcf_path), f'The input VCF file "{args.in_vcf_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.out_vcf_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'
    assert args.vep_script == 'vep' or os.path.isfile(args.vep_script), \
        f'The VEP script "{args.vep_script}" is not a binary or an invalid path.'
    assert 1 <= args.num_proc <= mp.cpu_count(), \
        f'Invalid number of processes "{args.num_proc:,d}". It must be in the range [1, {mp.cpu_count()}].'


def split_vcf_by_chrom(vcf_file_path: str) -> list:
    """ Split the input VCF file by each chromosome and return a list of split VCF file paths
    Each split file has a filename that ends with .{chromosome ID}.vcf.

    :param vcf_file_path: Input VCF file path
    :return: List of file paths of split VCFs
    """
    chrom_to_file = {}
    vcf_file_paths = []

    with open(vcf_file_path) as in_vcf_file:
        header = in_vcf_file.readline()

        for line in in_vcf_file:
            chrom = line.split('\t')[0]
            chr_vcf_file = chrom_to_file.get(chrom)

            if chr_vcf_file is None:
                chr_vcf_file_path = vcf_file_path.replace('.vcf', f'.{chrom}.vcf')
                vcf_file_paths.append(chr_vcf_file_path)
                chr_vcf_file = open(chr_vcf_file_path, 'w')
                chr_vcf_file.write(header)
                chrom_to_file[chrom] = chr_vcf_file

            chr_vcf_file.write(line)

    for chr_vcf_file in chrom_to_file.values():
        chr_vcf_file.close()

    return vcf_file_paths


def make_vep_cmd(vep_script: str, in_vcf_path: str, out_vcf_path: str, custom_path_dict: dict = None) -> str:
    """ Make a command to execute VEP and return it.
    """
    # Basic information
    cmd_args = [
        vep_script,
        '--assembly', 'GRCh38',
        '--offline',
        '--force_overwrite',
        '--format', 'vcf',
        '-i', in_vcf_path,
        '-o', out_vcf_path,
        '--vcf',
        '--no_stats',
        '--polyphen p',
    ]

    # Only the most severe consequence per gene.
    cmd_args += [
        '--per_gene',
        '--pick',
        '--pick_order', 'canonical,appris,tsl,biotype,ccds,rank,length',
    ]

    # Add options for nearest and distance
    cmd_args += [
        '--distance', '2000',
        '--nearest', 'symbol',
        '--symbol',
    ]

    # Add custom annotations
    if custom_path_dict is not None:
        for custom_file_key in custom_path_dict:
            custom_file_path = custom_path_dict[custom_file_key]

            if custom_file_path.endswith('vcf') or custom_file_path.endswith('vcf.gz'):
                cmd_args += [
                    '--custom',
                    ','.join([custom_file_path, custom_file_key, 'vcf', 'exact', '0', 'AF']),
                ]
            elif custom_file_path.endswith('bed') or custom_file_path.endswith('bed.gz'):
                cmd_args += [
                    '--custom',
                    ','.join([custom_file_path, custom_file_key, 'bed', 'overlap', '0']),
                ]
            elif custom_file_path.endswith('bw'):
                cmd_args += [
                    '--custom',
                    ','.join([custom_file_path, custom_file_key, 'bigwig', 'overlap', '0']),
                ]

    cmd = ' '.join(cmd_args)
    return cmd


def concat_vcf_files(out_vcf_path: str, vcf_file_paths: list):
    """ Concatenate the input VCF files into one VCF file

    :param out_vcf_path: Path of out VCF file
    :param vcf_file_paths: Paths of input VCF files
    """
    with open(out_vcf_path, 'w') as outfile:
        for i, vcf_file_path in enumerate(vcf_file_paths):
            with open(vcf_file_path) as vcf_file:
                if i == 0:  # Write all lines from the file
                    for line in vcf_file:
                        outfile.write(line)
                else:  # Write all lines except headers
                    for line in vcf_file:
                        if not line.startswith('#'):
                            outfile.write(line)


def annotate_by_bed(in_vcf_gz_path: str, out_vcf_path: str, annot_bed_path: str):
    """ Annotate variants in the input VCF file using the prepared annotation BED file """
    chroms = [f'chr{n}' for n in range(1, 23)]

    with pysam.TabixFile(in_vcf_gz_path) as in_vcf_file, pysam.TabixFile(annot_bed_path) as annot_bed_file, \
            open(out_vcf_path, 'w') as out_vcf_file:
        # Make and write headers
        vcf_headers = in_vcf_file.header
        annot_key_str = annot_bed_file.header[0].split('=')[1]
        annot_info_header = f'##INFO=<ID=ANNOT,Key={annot_key_str}>'
        vcf_headers.append(annot_info_header)
        vcf_headers[-1], vcf_headers[-2] = vcf_headers[-2], vcf_headers[-1]  # Swap

        for vcf_header in vcf_headers:
            print(vcf_header, file=out_vcf_file)

        # Annotate by the input BED file
        for chrom in chroms:
            var_iter = in_vcf_file.fetch(chrom, parser=pysam.asTuple())
            bed_iter = annot_bed_file.fetch(chrom, parser=pysam.asTuple())
            bed_memory = deque()
            variant = next(var_iter, None)

            while variant is not None:
                var_pos = int(variant[1]) - 1  # 1-based -> 0-based
                var_ref = variant[3]
                var_alt = variant[4]

                # Determine a search region for BED coordinates
                if len(var_ref) == 1:  # Insertion or substitution
                    region_start = var_pos
                    region_end = var_pos + 2 if len(var_alt) > 1 else var_pos + 1
                else:  # Deletion
                    region_start = var_pos + 1
                    region_end = region_start + len(var_ref) - 1

                # Get an annotation integer
                annot_int = 0
                stop_bed_iter = False

                # 1. Check the memory of previously checked BED coordinates
                while len(bed_memory) > 0 and int(bed_memory[0][2]) <= region_start:
                    # Remove non-overlapped coordinates
                    bed_memory.popleft()

                for bed in bed_memory:
                    if int(bed[1]) < region_end:  # Overlap
                        annot_int |= int(bed[3])
                    else:
                        stop_bed_iter = True
                        break

                # 2. Continuously iterate over the list BED coordinates and check
                if not stop_bed_iter:
                    bed = next(bed_iter, None)

                    while bed is not None:
                        bed_start = int(bed[1])
                        bed_end = int(bed[2])

                        if region_start < bed_end:
                            bed_memory.append(bed)
                            if bed_start < region_end:  # Overlap
                                annot_int |= int(bed[3])
                            else:
                                break

                        bed = next(bed_iter, None)

                print(str(variant) + f';ANNOT={annot_int}', file=out_vcf_file)
                variant = next(var_iter, None)


if __name__ == "__main__":
    main()

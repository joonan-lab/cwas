#!/usr/bin/env python
"""
Script for genomic and functional annotations using Variant Effect Predictor (VEP).
"""

import argparse
import multiprocessing as mp
import os

import yaml


def main():
    # Paths to essential configuration files
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    vep_custom_data_dir = os.path.join(project_dir, 'data', 'vep')
    vep_custom_conf_path = os.path.join(project_dir, 'conf', 'vep_custom_annotations.yaml')
    custom_file_path_dict = {}

    # Parse the configuration file
    with open(vep_custom_conf_path) as vep_custom_conf:
        custom_filename_dict = yaml.safe_load(vep_custom_conf)

        for custom_annot in custom_filename_dict:
            custom_file_path_dict[custom_annot] = os.path.join(vep_custom_data_dir, custom_filename_dict[custom_annot])

    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()
    print_args(args)
    check_args_validity(args)
    print()

    # Split the input file for each single chromosome
    if args.split_vcf:
        chr_vcf_file_paths = split_vcf_by_chrom(args.in_vcf_path)
        chr_vep_vcf_paths = []
        cmds = []

        # Make commands for VEP
        for chr_vcf_file_path in chr_vcf_file_paths:
            chr_vep_vcf_path = chr_vcf_file_path.replace('.vcf', '.vep.vcf')
            chr_vep_vcf_paths.append(chr_vep_vcf_path)
            cmd = make_vep_cmd(chr_vcf_file_path, chr_vep_vcf_path, custom_file_path_dict)
            cmds.append(cmd)

        # Run VEP in parallel
        pool = mp.Pool(args.num_proc)
        pool.map(os.system, cmds)
        pool.close()
        pool.join()

        # Concatenate the VEP outputs into one
        concat_vcf_files(args.out_vcf_path, chr_vep_vcf_paths)

        # Delete the split VCF files
        for chr_vcf_file_path in chr_vcf_file_paths:
            os.remove(chr_vcf_file_path)

        for chr_vep_vcf_path in chr_vep_vcf_paths:
            os.remove(chr_vep_vcf_path)
    else:
        cmd = make_vep_cmd(args.in_vcf_path, args.out_vcf_path, custom_file_path_dict)
        os.system(cmd)


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True, type=str, help='Input VCF file')
    parser.add_argument('-o', '--outfile', dest='out_vcf_path', required=False, type=str,
                        help='Path of the VCF output', default='vep_output.vcf')
    parser.add_argument('-s', '--split', dest='split_vcf', required=False, type=int, choices={0, 1},
                        help='Split the input VCF by chromosome and run VEP for each split VCF (Default: 0)',
                        default=0)
    parser.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                        help='Number of processes for this script (only necessary for split VCF files) (Default: 1)',
                        default=1)

    return parser


def print_args(args: argparse.Namespace):
    print(f'[Setting] The input VCF file: {args.in_vcf_path}')
    print(f'[Setting] The output path: {args.out_vcf_path}')

    if args.split_vcf:
        print(f'[Setting] Split the input VCF file by chromosome and run vep for each split VCF')
        print(f'[Setting] No. processes for this script: {args.num_proc:,d}')


def check_args_validity(args: argparse.Namespace):
    assert os.path.isfile(args.in_vcf_path), f'The input VCF file "{args.in_vcf_path}" cannot be found.'
    outfile_dir = os.path.dirname(args.out_vcf_path)
    assert outfile_dir == '' or os.path.isdir(outfile_dir), f'The outfile directory "{outfile_dir}" cannot be found.'
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


def make_vep_cmd(in_vcf_path: str, out_vcf_path: str, custom_file_path_dict: dict) -> str:
    """ Make a command to execute VEP and return it.
    """
    # Basic information
    cmd_args = [
        'vep',
        '--assembly', 'GRCh38',
        '--cache',
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
    for custom_annot in custom_file_path_dict:
        custom_file_path = custom_file_path_dict[custom_annot]

        if custom_file_path.endswith('vcf') or custom_file_path.endswith('vcf.gz'):
            cmd_args += [
                '--custom',
                ','.join([custom_file_path, custom_annot, 'vcf', 'exact', '0', 'AF']),
            ]
        elif custom_file_path.endswith('bw') or custom_file_path.endswith('bw.gz'):
            cmd_args += [
                '--custom',
                ','.join([custom_file_path, custom_annot, 'bigwig', 'overlap', '0']),
            ]
        else:  # BED files (.bed or .bed.gz)
            cmd_args += [
                '--custom',
                ','.join([custom_file_path, custom_annot, 'bed', 'overlap', '0']),
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


if __name__ == "__main__":
    main()

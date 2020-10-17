#!/usr/bin/env python
"""
Script for multiprocessing category-wide assocation study (CWAS)
"""
import argparse
import multiprocessing as mp
import os
from glob import glob

import yaml

from utils import get_curr_time, execute_cmd


def main():
    print(__doc__)

    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()
    print_args(args)

    # Path settings
    project_dir = os.path.dirname(os.path.abspath('.'))
    cwas_path_conf_path = os.path.join(project_dir, 'conf', 'cwas_paths.yaml')

    with open(cwas_path_conf_path, 'r') as cwas_path_conf:
        cwas_path_dict = yaml.safe_load(cwas_path_conf)

    # Make CMDs
    cwas_script = cwas_path_dict[args.step]
    cmds = []

    if args.step == 'annotate':
        infile_paths = sorted(glob(f'{args.in_dir}/*.vcf'))
        infile_paths = infile_paths[args.start_idx:args.end_idx]
        os.makedirs(args.out_dir, exist_ok=True)
        outfile_paths = [f'{args.out_dir}/{os.path.basename(in_vcf_path).replace(".vcf", ".annot.vcf")}'
                         for in_vcf_path in infile_paths]

        for infile_path, outfile_path in zip(infile_paths, outfile_paths):
            if args.force_overwrite or not os.path.isfile(outfile_path):
                cmd = f'{cwas_script} -i {infile_path} -o {outfile_path} --vep {args.vep_script};'
                cmds.append(cmd)

    elif args.step == 'categorize':
        infile_paths = sorted(glob(f'{args.in_dir}/*.vcf'))
        infile_paths = infile_paths[args.start_idx:args.end_idx]
        os.makedirs(args.out_dir, exist_ok=True)
        outfile_paths = [f'{args.out_dir}/{os.path.basename(in_vcf_path).replace(".vcf", ".cat_result.txt")}'
                         for in_vcf_path in infile_paths]

        for infile_path, outfile_path in zip(infile_paths, outfile_paths):
            if args.force_overwrite or not os.path.isfile(outfile_path):
                cmd = f'{cwas_script} -i {infile_path} -o {outfile_path};'
                cmds.append(cmd)

    elif args.step == 'burden_test':
        infile_paths = sorted(glob(f'{args.in_dir}/*.txt'))
        infile_paths = infile_paths[args.start_idx:args.end_idx]
        os.makedirs(args.out_dir, exist_ok=True)
        outfile_paths = [f'{args.out_dir}/{os.path.basename(in_vcf_path).replace(".txt", ".burden.txt")}'
                         for in_vcf_path in infile_paths]

        for infile_path, outfile_path in zip(infile_paths, outfile_paths):
            if args.force_overwrite or not os.path.isfile(outfile_path):
                cmd = f'{cwas_script} binom -i {infile_path} -o {outfile_path} -s {args.sample_file_path}'
                cmd += f' -a {args.adj_file_path};' if args.adj_file_path else ';'
                cmds.append(cmd)

    # Execute
    if args.num_proc == 1:
        for cmd in cmds:
            execute_cmd(cmd)
    else:
        pool = mp.Pool(args.num_proc)
        pool.map(execute_cmd, cmds)
        pool.close()
        pool.join()

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    parser = argparse.ArgumentParser(description=__doc__)

    def add_common_args(subparser: argparse.ArgumentParser):
        """ Add common arguments to the subparser """
        subparser.add_argument(
            '-i', '--in_dir', dest='in_dir', required=True, type=str, help='Directory of input files'
        )
        subparser.add_argument(
            '-o', '--out_dir', dest='out_dir', required=True, type=str, help='Directory of output VCFs'
        )
        subparser.add_argument(
            '-p', '--num_proc', dest='num_proc', required=False, type=int,
            help='Number of processes for this script (Default: 1)', default=1
        )
        subparser.add_argument(
            '--start_idx', dest='start_idx', required=False, type=int, default=0,
            help='Start index of a list of file paths (0-based) (Default: 0)'
        )
        subparser.add_argument(
            '--end_idx', dest='end_idx', required=False, type=int,
            help='End index of a list of file paths (Default: None)'
        )
        subparser.add_argument(
            '-f', '--force_overwrite', dest='force_overwrite', action='store_const',
            const=1, default=0, help='Force to overwrite when downloading data (Default: 0)'
        )

    subparsers = parser.add_subparsers(description='A name of a step of CWAS {annotate, categorize, burden_test}',
                                       dest='step')
    parser_annot = subparsers.add_parser(
        'annotate',
        description='Multiprocessing variant annotation in CWAS',
        help='Multiprocessing variant annotation in CWAS (arg "annotate -h" for usage)'
    )
    parser_annot.add_argument('--vep', dest='vep_script', required=False, type=str,
                              help='Path of a Perl script to execute VEP (Default: vep (binary))', default='vep')
    add_common_args(parser_annot)

    parser_cat = subparsers.add_parser(
        'categorize',
        description='Multiprocessing variant annotation in CWAS',
        help='Multiprocessing variant categorization in CWAS (arg "categorize -h" for usage)'
    )
    add_common_args(parser_cat)

    parser_burden = subparsers.add_parser(
        'burden_test',
        description='Multiprocessing burden binomial tests in CWAS',
        help='Multiprocessing burden binomial tests in CWAS (arg "burden_test -h" for usage)'
    )
    parser_burden.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                               help='File listing sample IDs with their families and sample_types (case or ctrl)')
    parser_burden.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                               help='File that contains adjustment factors for No. DNVs of each sample', default='')
    add_common_args(parser_burden)

    return parser


def print_args(args: argparse.Namespace):
    print(f'[Setting] CWAS step: {args.step}')
    print(f'[Setting] Input directory: {args.in_dir}')
    print(f'[Setting] Output directory: {args.out_dir}')
    print(f'[Setting] No. processes for this script: {args.num_proc:,d}')
    print(f'[Setting] Start index of file paths: {args.start_idx}')
    print(f'[Setting] End index of file paths: {args.end_idx}')

    if args.step == 'annotate':
        print(f'[Setting] VEP script: {args.vep_script}')


if __name__ == '__main__':
    main()

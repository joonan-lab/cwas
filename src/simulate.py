#!/usr/bin/env python
"""
Script for simulating the generation of random mutations.
The random mutations will be used to get simulated p-values.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import os

import yaml

from utils import get_curr_time


def main():
    # Paths for this script
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    tmp_dir = os.path.join(project_dir, 'tmp')
    filepath_conf_path = os.path.join(project_dir, 'conf', 'filepaths.yaml')
    fileurl_conf_path = os.path.join(project_dir, 'conf', 'fileurls.yaml')
    os.makedirs(tmp_dir, exist_ok=True)

    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()

    if args.mode == 'download':
        # Download essential data
        with open(filepath_conf_path) as filepath_conf_file:
            filepath_conf = yaml.safe_load(filepath_conf_file)
            filepath_dict = filepath_conf['simulate']

        with open(fileurl_conf_path) as fileurl_conf_file:
            fileurl_conf = yaml.safe_load(fileurl_conf_file)
            fileurl_dict = fileurl_conf['simulate']

        data_dir = os.path.join(project_dir, filepath_dict['data_dir'])
        os.makedirs(data_dir, exist_ok=True)

        for data_key in fileurl_dict:
            data_dest_path = os.path.join(project_dir, filepath_dict[data_key])

            if os.path.isfile(data_dest_path):
                print(f'[INFO] "{data_dest_path}" already exists.')
            else:
                cmd = f'wget -O {data_dest_path} {fileurl_dict[data_key]}'
                print(f'[CMD] {cmd}')
                os.system(cmd)


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    # Create a top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(description='Execution mode', dest='mode')

    # Create a parser for downloading data
    subparsers.add_parser('download', description='Download data essential for the simulation',
                          help='Download data (arg "download -h" for usage)')

    # Create a parser for processing the downloaded data as a preparation step
    subparsers.add_parser('prepare', description='Process the downloaded data for preparation',
                          help='Process the downloaded data (arg "prepare -h" for usage')

    # Create a parser for generating random mutations
    parser_mut = subparsers.add_parser('mutation', description='Generate random mutations',
                                       help='Generate random mutations (arg "mutation -h" for usage')
    parser_mut.add_argument('-i', '--infile', dest='in_vcf_path', required=True, type=str,
                            help='Input VCF file which is referred to generate random mutations')
    parser_mut.add_argument('-o', '--out_dir', dest='out_dir', required=True, type=str,
                            help='Directory of outputs that lists random mutations. '
                                 'The number of outputs will be the same with the number of simulations.')
    parser_mut.add_argument('-n', '--num_sim', dest='num_sim', required=False, type=int,
                            help='Number of simulations to generate random mutations (Default: 1)', default=1)

    return parser


if __name__ == '__main__':
    main()

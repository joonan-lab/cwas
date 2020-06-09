#!/usr/bin/env python
"""
Script for simulating the generation of random mutations.
The random mutations will be used to get simulated p-values.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import gzip
import os

import yaml

from utils import get_curr_time
from collections import defaultdict


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

    if args.mode == 'download':  # Download essential data
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

    elif args.mode == 'prepare':  # Process the downloaded data
        with open(filepath_conf_path) as filepath_conf_file:
            filepath_conf = yaml.safe_load(filepath_conf_file)
            filepath_dict = filepath_conf['simulate']

        # path settings
        gap_path = os.path.join(project_dir, filepath_dict['gap'])
        lcr_path = os.path.join(project_dir, filepath_dict['lcr'])
        sort_gap_path = os.path.join(project_dir, filepath_dict['sort_gap'])
        mask_region_path = os.path.join(project_dir, filepath_dict['mask_region'])
        chrom_size_path = os.path.join(project_dir, filepath_dict['chrom_size'])

        if not os.path.isfile(sort_gap_path):
            cmd = f'gunzip -c {gap_path} | cut -f2,3,4 - | sort -k1,1 -k2,2n | gzip > {sort_gap_path};'
            print(f'[{get_curr_time()}, Progress] Sort the gap regions')
            os.system(cmd)

        if not os.path.isfile(mask_region_path):
            cmd = f'zcat {sort_gap_path} {lcr_path} | sortBed -i stdin | gzip > {mask_region_path}'
            print(f'[{get_curr_time()}, Progress] Merge the gap and LCR regions')
            os.system(cmd)

        chroms = [f'chr{n}' for n in range(1, 23)]
        for chrom in chroms:
            in_fa_gz_path = os.path.join(project_dir, filepath_dict[chrom])
            in_fa_path = in_fa_gz_path.replace('.gz', '')
            out_fa_path = os.path.join(project_dir, filepath_dict[f'{chrom}_masked'])

            if not os.path.isfile(out_fa_path):
                print(f'[{get_curr_time()}, Progress] Mask the {chrom} fasta file and index the output')
                cmd = f'gunzip {in_fa_gz_path};'
                cmd += f'maskFastaFromBed -fi {in_fa_path} -fo {out_fa_path} -bed {mask_region_path};'
                cmd += f'samtools faidx {out_fa_path};'
                os.system(cmd)

        if not os.path.isfile(chrom_size_path):
            print(f'[{get_curr_time()}, Progress] Calculate mapped, AT/GC, and effective sizes of each chromosome')

            with open(chrom_size_path, 'w') as outfile:
                print('Chrom', 'Size', 'Mapped', 'AT', 'GC', 'Effective', sep='\t', file=outfile)

                for chrom in chroms:
                    print(f'[{get_curr_time()}, Progress] {chrom}')
                    mask_fa_path = os.path.join(project_dir, filepath_dict[f'{chrom}_masked'])
                    fa_idx_path = mask_fa_path.replace('.gz', '.fai')

                    with open(fa_idx_path) as fa_idx_file:
                        line = fa_idx_file.read()
                        fields = line.split('\t')
                        chrom_size = int(fields[1])
                        base_start_idx = int(fields[2])

                    with gzip.open(mask_fa_path, 'rt') as mask_fa_file:
                        base_cnt_dict = defaultdict(int)
                        mask_fa_file.seek(base_start_idx)
                        seq = mask_fa_file.read()

                        for base in seq:
                            base_cnt_dict[base.upper()] += 1

                    map_size = chrom_size - base_cnt_dict['N']
                    at_size = base_cnt_dict['A'] + base_cnt_dict['T']
                    gc_size = base_cnt_dict['G'] + base_cnt_dict['C']
                    effect_size = (at_size / 1.75) + gc_size
                    print(chrom, chrom_size, map_size, at_size, gc_size, effect_size, sep='\t', file=outfile)

        print(f'[{get_curr_time()}, Progress] Done')

    elif args.mode == 'mutation':  # Generate random mutations
        pass

    else:
        pass


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
    parser_mut.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                            help='File listing sample IDs with their families and sample_types (case or ctrl)')
    parser_mut.add_argument('-o', '--out_dir', dest='out_dir', required=True, type=str,
                            help='Directory of outputs that lists random mutations. '
                                 'The number of outputs will be the same with the number of simulations.')
    parser_mut.add_argument('-n', '--num_sim', dest='num_sim', required=False, type=int,
                            help='Number of simulations to generate random mutations (Default: 1)', default=1)

    return parser


if __name__ == '__main__':
    main()

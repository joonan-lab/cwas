#!/usr/bin/env python
"""
Script for simulating the generation of random mutations.
The random mutations will be used to get simulated p-values.

For more detailed information, please refer to An et al., 2018 (PMID 30545852).

"""
import argparse
import gzip
import multiprocessing as mp
import os
from collections import defaultdict
from functools import partial

import numpy as np
import pandas as pd
import yaml

from fastafile import FastaFile
from utils import get_curr_time, div_list


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
        print(f'[{get_curr_time()}, Progress] Load and parse the data to generate random mutations')
        # Load the input data
        variant_df = parse_vcf(args.in_vcf_path)
        sample_df = pd.read_table(args.sample_file_path, index_col='SAMPLE')

        # Load data from the preparation step
        with open(filepath_conf_path) as filepath_conf_file:
            filepath_conf = yaml.safe_load(filepath_conf_file)
            filepath_dict = filepath_conf['simulate']

        chrom_size_path = os.path.join(project_dir, filepath_dict['chrom_size'])

        if not os.path.isfile(chrom_size_path):
            raise FileNotFoundError(f'"{chrom_size_path}" does not exist. Run the "prepare" step first.')

        chrom_size_df = pd.read_table(chrom_size_path)

        # Extract values from the DataFrames
        samples = variant_df['SAMPLE'].values
        refs = variant_df['REF'].values
        alts = variant_df['ALT'].values
        variant_labels = np.vectorize(label_variant)(refs, alts)
        sample_to_fam = sample_df.to_dict()['FAMILY']
        chrom_eff_sizes = chrom_size_df['Effective'].values
        chrom_probs = chrom_eff_sizes / np.sum(chrom_eff_sizes)  # Normalization
        chrom_sizes = chrom_size_df['Size'].values
        chroms = chrom_size_df['Chrom'].values

        # Make dictionaries to generate random mutation
        fam_to_label_cnt = {}  # Key: Family ID, Value: Array of label counts (Each index corresponds to each label.)
        fam_to_sample_set = {}  # Key: Family ID, Value: Set of sample IDs available in the input VCF file

        for sample, variant_label in zip(samples, variant_labels):
            family = sample_to_fam[sample]
            label_cnt_arr = fam_to_label_cnt.get(family, np.zeros(4, dtype=int))
            label_cnt_arr[variant_label] += 1
            fam_to_label_cnt[family] = label_cnt_arr

            sample_set = fam_to_sample_set.get(family, set())
            sample_set.add(sample)
            fam_to_sample_set[family] = sample_set

        # Make a dictionary for FASTA file paths masked in the previous preparation step.
        unq_chroms = np.unique(chroms)
        fasta_path_dict = {}

        for chrom in unq_chroms:
            fasta_file_path = os.path.join(project_dir, filepath_dict[f'{chrom}_masked'])
            fasta_path_dict[chrom] = fasta_file_path

        # Make files listing random mutations
        print(f'[{get_curr_time()}, Progress] Make files listing generated random mutations')
        os.makedirs(args.out_dir, exist_ok=True)
        output_paths = []

        for n in range(args.num_sim):
            output_filename = f'{args.out_tag}.{n + 1:05d}.vcf'
            output_path = os.path.join(args.out_dir, output_filename)
            output_paths.append(output_path)

        if args.num_proc == 1:
            make_rand_mut_files(output_paths, fam_to_label_cnt, fam_to_sample_set, fasta_path_dict,
                                chrom_probs, chrom_sizes)
        else:
            output_paths_subs = div_list(output_paths, args.num_proc)
            pool = mp.Pool(args.num_proc)
            pool.map(
                partial(make_rand_mut_files,
                        fam_to_label_cnt=fam_to_label_cnt,
                        fam_to_sample_set=fam_to_sample_set,
                        fasta_path_dict=fasta_path_dict,
                        chrom_probs=chrom_probs,
                        chrom_sizes=chrom_sizes,
                        ),
                output_paths_subs
            )
            pool.close()
            pool.join()

        print(f'[{get_curr_time()}, Progress] Done')


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
    parser_mut.add_argument('-o', '--out_dir', dest='out_dir', required=False, type=str,
                            help='Directory of outputs that lists random mutations. '
                                 'The number of outputs will be the same with the number of simulations. '
                                 '(Default: ./random-mutation)', default='random-mutation')
    parser_mut.add_argument('-t', '--out_tag', dest='out_tag', required=False, type=str,
                            help='Prefix of output files. Each output file name will start with this tag. '
                                 '(Default: rand_mut)', default='rand_mut')
    parser_mut.add_argument('-n', '--num_sim', dest='num_sim', required=False, type=int,
                            help='Number of simulations to generate random mutations (Default: 1)', default=1)
    parser_mut.add_argument('-p', '--num_proc', dest='num_proc', required=False, type=int,
                            help='Number of processes for this script (only necessary for split VCF files) '
                                 '(Default: 1)', default=1)

    return parser


def parse_vcf(vcf_path: str, rdd_colnames: list = None) -> pd.DataFrame:
    """ Parse the VCF file and make a pandas.DataFrame object listing the annotated variants.

    :param vcf_path: The path of the VCF file
    :param rdd_colnames: The list of column names redundant for CWAS
                         (Warning: Unavailable column names will be ignored.)
    :return: The DataFrame object listing annotated variants
    """
    variant_df_rows = []
    variant_df_colnames = []

    # Parse the VCF file
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):  # The comments
                if line.startswith('#CHROM'):  # The header
                    variant_df_colnames = line[1:].rstrip('\n').split('\t')
            else:
                variant_df_row = line.rstrip('\n').split('\t')
                variant_df_rows.append(variant_df_row)

    vcf_df = pd.DataFrame(variant_df_rows, columns=variant_df_colnames)

    # Parse the INFO field
    info_strs = vcf_df['INFO'].values
    info_dicts = list(map(parse_info_str, info_strs))
    info_df = pd.DataFrame(info_dicts)

    # Concatenate those DataFrames
    variant_df = pd.concat([vcf_df.drop(columns='INFO'), info_df], axis='columns')

    # Trim the columns redundant for CWAS
    if rdd_colnames is not None:
        variant_df.drop(columns=rdd_colnames, inplace=True, errors='ignore')

    return variant_df


def parse_info_str(info_str: str) -> dict:
    """ Parse the string in the INFO field of the VCF file and make a dictionary """
    info_dict = {}
    key_value_pairs = info_str.split(';')

    for key_value_pair in key_value_pairs:
        key, value = key_value_pair.split('=', 1)
        info_dict[key] = value

    return info_dict


def label_variant(ref: str, alt: str) -> int:
    """ Return an integer according to the type of the input small variant
    """
    assert len(ref) == 1 or len(alt) == 1, 'This function current support only small variants such as SNV and INDEL.'

    if len(ref) == 1 and len(alt) == 1:
        return 0  # SNV
    elif len(ref) % 3 == 1 or len(alt) % 3 == 1:
        return 1  # Indel0
    elif len(ref) % 3 == 2 or len(alt) % 3 == 2:
        return 2  # Indel1
    else:  # len(ref) % 3 == 0 or len(alt) % 3 == 0:
        return 3  # Indel2


def make_rand_mut_files(output_paths: list, fam_to_label_cnt: dict, fam_to_sample_set: dict, fasta_path_dict: dict,
                        chrom_probs: np.ndarray, chrom_sizes: np.ndarray):
    """ Make VCF files listing random mutations
    This is a wrapper for the 'make_rand_mut_file' function to support multiprocessing.
    """
    # Open the FASTA files
    fasta_file_dict = {chrom: FastaFile(fasta_path_dict[chrom]) for chrom in fasta_path_dict}

    for output_path in output_paths:
        make_rand_mut_file(output_path, fam_to_label_cnt, fam_to_sample_set, fasta_file_dict, chrom_probs, chrom_sizes)

    # Close the FASTA files
    for chrom in fasta_file_dict:
        fasta_file_dict[chrom].close()


def make_rand_mut_file(output_path: str, fam_to_label_cnt: dict, fam_to_sample_set: dict, fasta_file_dict: dict,
                       chrom_probs: np.ndarray, chrom_sizes: np.ndarray):
    """ Make a VCF file listing random mutations """
    rand_variants = []

    for fam in fam_to_label_cnt:
        label_cnt_arr = fam_to_label_cnt[fam]
        sample_ids = list(fam_to_sample_set[fam])

        for label, label_cnt in enumerate(label_cnt_arr):
            for i in range(label_cnt):
                rand_variant = make_random_mutation(label, sample_ids, fasta_file_dict, chrom_probs, chrom_sizes)
                rand_variants.append(rand_variant)

    rand_variants.sort(key=lambda x: (x.get('chrom'), x.get('pos')))
    write_variant_list(output_path, rand_variants)


def make_random_mutation(label: int, sample_ids: list, fasta_file_dict: dict,
                         chrom_probs: np.ndarray, chrom_sizes: np.ndarray) -> dict:
    """ Generate and return a random mutation in the VCF format """
    sample_id = np.random.choice(sample_ids)
    ref = None
    alt = None

    while True:
        chrom_idx = np.random.choice(range(len(chrom_probs)), p=chrom_probs)
        chrom_size = chrom_sizes[chrom_idx]
        chrom = f'chr{chrom_idx + 1}'
        fasta_file = fasta_file_dict[chrom]

        pos = np.random.randint(chrom_size)
        base = fasta_file.get_base(chrom, pos).upper()

        if base == 'N':
            continue

        ref, alt = pick_mutation()

        if base != ref:
            continue

        break

    alt += 'A' * label

    variant = {
        'chrom': chrom,
        'pos': pos + 1,  # 0-based -> 1-based
        'id': f'{chrom}:{pos + 1}:{ref}:{alt}',
        'ref': ref,
        'alt': alt,
        'qual': '.',
        'filter': '.',
        'info': f'SAMPLE={sample_id}'
    }
    return variant


def pick_mutation() -> (str, str):
    """ Get a mutation from the mutation distribution model (Lynch 2010). """
    x = np.random.uniform(0, 1)

    if x <= 0.211062887:
        ref = 'C'
        alt = 'T'
    elif x <= 0.422125774:
        ref = 'G'
        alt = 'A'
    elif x <= 0.551200326:
        ref = 'A'
        alt = 'G'
    elif x <= 0.680274877:
        ref = 'T'
        alt = 'C'
    elif x <= 0.728393387:
        ref = 'G'
        alt = 'T'
    elif x <= 0.776511898:
        ref = 'C'
        alt = 'A'
    elif x <= 0.821985623:
        ref = 'G'
        alt = 'C'
    elif x <= 0.867459349:
        ref = 'C'
        alt = 'G'
    elif x <= 0.900590744:
        ref = 'T'
        alt = 'A'
    elif x <= 0.933722139:
        ref = 'A'
        alt = 'T'
    elif x <= 0.96686107:
        ref = 'A'
        alt = 'C'
    else:
        ref = 'T'
        alt = 'G'

    return ref, alt


def write_variant_list(out_vcf_path, variants):
    with open(out_vcf_path, 'w') as outfile:
        print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', sep='\t', file=outfile)

        for variant in variants:
            print(variant['chrom'], variant['pos'], variant['id'], variant['ref'], variant['alt'],
                  variant['qual'], variant['filter'], variant['info'], sep='\t', file=outfile)


if __name__ == '__main__':
    main()

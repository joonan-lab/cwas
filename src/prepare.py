#!/usr/bin/env python
"""
Script for data preparation for category-wide association study (CWAS)
#TODO: Multiprocessing, File existence check
"""
import argparse
import gzip
import os
from collections import defaultdict

import pyBigWig as pbw
import pysam
import yaml

from utils import get_curr_time


def main():
    # Paths for this script
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    ori_file_conf_path = os.path.join(project_dir, 'conf', 'download_filepaths.yaml')  # Files downloaded already
    target_file_conf_path = os.path.join(project_dir, 'conf', 'prepare_filepaths.yaml')  # Files that will be made

    # Parse the configuration files
    with open(ori_file_conf_path) as ori_file_conf_file, open(target_file_conf_path) as target_file_conf_file:
        ori_filepath_dict = yaml.safe_load(ori_file_conf_file)
        target_filepath_dict = yaml.safe_load(target_file_conf_file)

        for file_group in target_filepath_dict:
            for file_key in target_filepath_dict[file_group]:
                ori_filepath_dict[file_group][file_group] = \
                    os.path.join(project_dir, ori_filepath_dict[file_group][file_key])
                target_filepath_dict[file_group][file_group] = \
                    os.path.join(project_dir, target_filepath_dict[file_group][file_key])

    # Parse arguments
    parser = create_arg_parser()
    args = parser.parse_args()

    if args.step == 'simulate':
        print(f'[{get_curr_time()}, Progress] Prepare data for random mutation simulation')
        ori_filepath_dict = ori_filepath_dict['simulate']
        target_filepath_dict = target_filepath_dict['simulate']
        make_mask_region_bed(ori_filepath_dict, target_filepath_dict)
        mask_fasta(ori_filepath_dict, target_filepath_dict)
        make_chrom_size_txt(target_filepath_dict)

    elif args.step == 'annotate':
        print(f'[{get_curr_time()}, Progress] Prepare data for variant annotation')
        ori_filepath_dict = ori_filepath_dict['annotate']
        target_filepath_dict = target_filepath_dict['annotate']

        # Filter entries of Yale PsychENCODE BED files
        file_keys = ['Yale_H3K27ac_CBC', 'Yale_H3K27ac_DFC']
        for file_key in file_keys:
            filt_yale_bed(ori_filepath_dict[file_key], target_filepath_dict[file_key])
            bgzip_tabix(target_filepath_dict[file_key])

        # Make BED files for conservation scores from the BigWig files
        chrom_size_path = ori_filepath_dict['chrom_size']
        chrom_size_dict = {}

        with open(chrom_size_path) as chrom_size_file:
            for line in chrom_size_file:
                fields = line.strip().split('\t')
                chrom_size_dict[fields[0]] = int(fields[1])

        make_bed_from_bw(ori_filepath_dict['phyloP46wayVt'], target_filepath_dict['phyloP46wayVt'],
                         2.0, chrom_size_dict)
        make_bed_from_bw(ori_filepath_dict['phastCons46wayVt'], target_filepath_dict['phastCons46wayVt'],
                         0.2, chrom_size_dict)
        bgzip_tabix(target_filepath_dict['phyloP46wayVt'])
        bgzip_tabix(target_filepath_dict['phastCons46wayVt'])

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    # Create a top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(description='Step for which data is prepared', dest='step')
    subparsers.add_parser('simulate', description='Prepare data to simulate random mutations',
                          help='Prepare data to simulate random mutations (arg "simulate -h" for usage)')
    return parser


def make_mask_region_bed(ori_filepath_dict: dict, target_filepath_dict: dict):
    """ Make a bed file listing masked regions by merging gap and LCR regions """
    gap_path = ori_filepath_dict['gap']
    lcr_path = ori_filepath_dict['lcr']
    sort_gap_path = target_filepath_dict['gap']
    mask_region_path = target_filepath_dict['mask_region']

    if os.path.isfile(sort_gap_path):
        print(f'[{get_curr_time()}, Progress] There is already a list of sorted gap regions so skip this step')
    else:
        cmd = f'gunzip -c {gap_path} | cut -f2,3,4 - | sort -k1,1 -k2,2n | gzip > {sort_gap_path};'
        print(f'[{get_curr_time()}, Progress] Sort the gap regions')
        print(f'[{get_curr_time()} CMD] {cmd}')
        os.system(cmd)

    if os.path.isfile(mask_region_path):
        print(f'[{get_curr_time()}, Progress] Masked regions have already made so skip this step')
    else:
        cmd = f'zcat {sort_gap_path} {lcr_path} | sortBed -i stdin | gzip > {mask_region_path}'
        print(f'[{get_curr_time()}, Progress] Make masked regions by merging the gap and LCR regions')
        print(f'[{get_curr_time()} CMD] {cmd}')
        os.system(cmd)


def mask_fasta(ori_filepath_dict: dict, target_filepath_dict: dict):
    """ Mask regions on the genome files (FASTA files) and save results as new FASTA files """
    mask_region_path = target_filepath_dict['mask_region']
    chroms = [f'chr{n}' for n in range(1, 23)]

    for chrom in chroms:
        in_fa_gz_path = ori_filepath_dict[chrom]
        in_fa_path = in_fa_gz_path.replace('.gz', '')
        out_fa_path = target_filepath_dict[chrom]

        if os.path.isfile(out_fa_path):
            print(f'[{get_curr_time()}, Progress] Masked fasta file for {chrom} already exists so skip this step')
        else:
            print(f'[{get_curr_time()}, Progress] Mask the {chrom} fasta file and index the output')
            cmd = f'gunzip {in_fa_gz_path};'
            cmd += f'maskFastaFromBed -fi {in_fa_path} -fo {out_fa_path} -bed {mask_region_path};'
            cmd += f'samtools faidx {out_fa_path};'
            print(f'[{get_curr_time()} CMD] {cmd}')
            os.system(cmd)


def make_chrom_size_txt(target_filepath_dict: dict):
    """ Make a txt file listing total size, mapped, AT/GC, and effective sizes of each chromosome """
    chrom_size_path = target_filepath_dict['chrom_size']
    chroms = [f'chr{n}' for n in range(1, 23)]

    if os.path.isfile(chrom_size_path):
        print(f'[{get_curr_time()}, Progress] A file listing total, mapped, AT/GC, and effective sizes already exists '
              f'so skip this step')
    else:
        print(f'[{get_curr_time()}, Progress] Make a file listing total, mapped, AT/GC, and effective sizes '
              f'of each chromosome')

        with open(chrom_size_path, 'w') as outfile:
            print('Chrom', 'Size', 'Mapped', 'AT', 'GC', 'Effective', sep='\t', file=outfile)

            for chrom in chroms:
                print(f'[{get_curr_time()}, Progress] {chrom}')
                mask_fa_path = target_filepath_dict[f'{chrom}']
                fa_idx_path = mask_fa_path + '.fai'

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


def filt_yale_bed(in_bed_path: str, out_bed_path: str):
    """ Filter entries of a BED file from Yale (PsychENCODE Consortium) and
    write a BED file listing the filtered entries.
    """
    with pysam.TabixFile(in_bed_path) as in_bed_file, open(out_bed_path, 'w') as out_bed_file:
        for bed_fields in in_bed_file.fetch(parser=pysam.asTuple()):
            bed_val = max([int(x.split('_')[0]) for x in bed_fields[3].split('&')])

            if bed_val > 1:
                print(*bed_fields, sep='\t', file=out_bed_file)


def bgzip_tabix(bed_path: str):
    """ Block compression (bgzip) and make an index (tabix) """
    cmd = f'bgzip {bed_path};'
    cmd += f'tabix {bed_path + ".gz"};'
    print(f'[{get_curr_time()}, CMD] {cmd}')
    exit_val = os.system(cmd)

    if exit_val != 0:
        print(f'[{get_curr_time()}, WARNING] This CMD is failed with this exit value {exit_val}.')


def make_bed_from_bw(in_bw_path: str, out_bed_path: str, cutoff: float, chrom_size_dict: dict):
    """ Make a BED file from a BigWig file.
    A BED entry covers a region of which all positions have data values more than or equal to the input cutoff
    """
    chroms = [f'chr{n}' for n in range(1, 23)]
    bin_size = 1000000  # Size of chromosomal bins

    with open(out_bed_path, 'w') as bed_file:
        for chrom in chroms:
            bed_entries = make_bed_entries(in_bw_path, chrom, chrom_size_dict[chrom], bin_size, cutoff)

            for bed_entry in bed_entries:
                print(*bed_entry, 1, sep='\t', file=bed_file)


def make_bed_entries(bw_path: str, chrom: str, chrom_size: int, chrom_bin_size: int, cutoff: float) -> list:
    """ Make BED entries from the input BigWig file """
    chrom_bins = make_bins(chrom_bin_size, chrom_size)
    bed_entries = []
    interval_stack = []

    with pbw.open(bw_path) as bw_file:
        for start, end in chrom_bins:
            intervals = bw_file.intervals(chrom, start, end)

            if not intervals:
                continue

            for interval in intervals:
                if intervals[2] < cutoff:
                    continue

                if not interval_stack or interval_stack[-1][1] == intervals[0]:  # Continuous interval
                    interval_stack.append(interval)
                else:
                    bed_entries.append((chrom, interval_stack[0][0], interval_stack[-1][1]))
                    interval_stack = [interval]

    # Make a bed entry using remain intervals in the stack
    bed_entries.append((chrom, interval_stack[0][0], interval_stack[-1][1]))
    return bed_entries


def make_bins(bin_size: int, total_size: int) -> list:
    bins = []
    bin_cnt = total_size // bin_size
    remain = total_size % bin_size

    for i in range(bin_cnt):
        bins.append((bin_size * i, bin_size * (i + 1)))

    if remain != 0:
        bins.append((bin_cnt * bin_size, bin_cnt * bin_size + remain))

    return bins


if __name__ == '__main__':
    main()

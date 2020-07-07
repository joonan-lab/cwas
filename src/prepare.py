#!/usr/bin/env python
"""
Script for data preparation for category-wide association study (CWAS)
#TODO: Multiprocessing
"""
import argparse
import gzip
import os
from collections import defaultdict

import numpy as np
import pyBigWig as pbw
import pysam
import yaml

from utils import get_curr_time, execute_cmd


def main():
    # Paths for this script
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    ori_file_conf_path = os.path.join(project_dir, 'conf', 'download_filepaths.yaml')  # Files downloaded already
    target_file_conf_path = os.path.join(project_dir, 'conf', 'prepare_filepaths.yaml')  # Files that will be made

    # Parse the configuration files
    ori_filepath_dict = {}
    target_filepath_dict = {}

    with open(ori_file_conf_path) as ori_file_conf_file, open(target_file_conf_path) as target_file_conf_file:
        ori_filepath_conf = yaml.safe_load(ori_file_conf_file)
        target_filepath_conf = yaml.safe_load(target_file_conf_file)

        for file_group in ori_filepath_conf:
            ori_filepath_dict[file_group] = {}

            for file_key in ori_filepath_conf[file_group]:
                ori_filepath_dict[file_group][file_key] = \
                    os.path.join(project_dir, ori_filepath_conf[file_group][file_key])

        for file_group in target_filepath_conf:
            target_filepath_dict[file_group] = {}

            for file_key in target_filepath_conf[file_group]:
                target_filepath_dict[file_group][file_key] = \
                    os.path.join(project_dir, target_filepath_conf[file_group][file_key])

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

        # Path settings for this step
        annot_dir = os.path.join(project_dir, 'data', 'annotate')
        annot_conf_path = os.path.join(project_dir, 'conf', 'annotation.yaml')
        annot_bed_path_dict = {}  # Dictionary of custom BED file paths

        # Parse the configuration file for this step
        with open(annot_conf_path) as annot_conf_file:
            annot_filename_dict = yaml.safe_load(annot_conf_file)

            for annot_key in annot_filename_dict:
                annot_bed_path_dict[annot_key] = os.path.join(annot_dir, annot_filename_dict[annot_key])

        # Merge annotation information of the custom BED files
        merge_bed_path = os.path.join(annot_dir, 'merged_annotation.bed')
        merge_annot(merge_bed_path, annot_bed_path_dict)

    print(f'[{get_curr_time()}, Progress] Done')


def create_arg_parser() -> argparse.ArgumentParser:
    """ Create an argument parser for this script and return it """
    # Create a top-level argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(description='Step for which data is prepared', dest='step')
    subparsers.add_parser('simulate', description='Prepare data to simulate random mutations',
                          help='Prepare data to simulate random mutations (arg "simulate -h" for usage)')
    subparsers.add_parser('annotate', description='Prepare data to annotate variants',
                          help='Prepare data to annotate variants (arg "annotate -h" for usage)')
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
        execute_cmd(cmd)

    if os.path.isfile(mask_region_path):
        print(f'[{get_curr_time()}, Progress] Masked regions have already made so skip this step')
    else:
        cmd = f'zcat {sort_gap_path} {lcr_path} | sortBed -i stdin | gzip > {mask_region_path}'
        print(f'[{get_curr_time()}, Progress] Make masked regions by merging the gap and LCR regions')
        execute_cmd(cmd)


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
            execute_cmd(cmd)


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
    if os.path.isfile(out_bed_path) or os.path.isfile(out_bed_path + '.gz'):
        print(f'[{get_curr_time()}, Progress] Filtered bed file "{out_bed_path}" already exists '
              f'so skip this filtering step')
    else:
        print(f'[{get_curr_time()}, Progress] Filter the bed file "{in_bed_path}"')
        with pysam.TabixFile(in_bed_path) as in_bed_file, open(out_bed_path, 'w') as out_bed_file:
            for bed_fields in in_bed_file.fetch(parser=pysam.asTuple()):
                bed_val = max([int(x.split('_')[0]) for x in bed_fields[3].split('&')])

                if bed_val > 1:
                    print(*bed_fields, sep='\t', file=out_bed_file)


def bgzip_tabix(bed_path: str):
    """ Block compression (bgzip) and make an index (tabix) """
    bed_gz_path = bed_path + '.gz'

    if os.path.isfile(bed_gz_path):
        print(f'[{get_curr_time()}, Progress] A bgzipped file for "{bed_path}" already exists so skip this step')
    else:
        print(f'[{get_curr_time()}, Progress] bgzip and tabix for "{bed_path}"')
        cmd = f'bgzip {bed_path};'
        cmd += f'tabix {bed_gz_path};'
        execute_cmd(cmd)


def make_bed_from_bw(in_bw_path: str, out_bed_path: str, cutoff: float, chrom_size_dict: dict):
    """ Make a BED file from a BigWig file.
    A BED entry covers a region of which all positions have data values more than or equal to the input cutoff
    """
    chroms = [f'chr{n}' for n in range(1, 23)]
    bin_size = 1000000  # Size of chromosomal bins

    if os.path.isfile(out_bed_path) or os.path.isfile(out_bed_path + '.gz'):
        print(f'[{get_curr_time()}, Progress] A BED file for "{in_bw_path}" already exists so skip this step')
    else:
        print(f'[{get_curr_time()}, Progress] Make a BED file for "{in_bw_path}"')
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


def merge_annot(out_bed_path: str, bed_path_dict: dict):
    """ Merge all annotation information of coordinates from all the annotation BED files into one file.

    For example, assume that following two coordinates exist in bed file A and B, respectively.
    - chr1  100 300 1
    - chr1  200 400 1

    If these coordinates are merged, than following coordinates will be newly produced.
    - chr1  100 200 annotation A
    - chr1  200 300 annotation A & B
    - chr1  300 400 annotation B

    The information of annotation for each coordinate will be represented as an annotation integer, which is a one-hot
    encoding of annotation information. For example, if all of the annotation files are A, B, and C, then an annotation
    integer of the second coordinate above is represented as 6(b'110).

    So here is a final result.
    - chr1  100 200 4
    - chr1  200 300 6
    - chr1  300 400 2
    """
    if os.path.isfile(out_bed_path):
        print(f'[{get_curr_time()}, Progress] A annotation-merged bed file already exists so skip the merge step.')
    else:
        print(f'[{get_curr_time()}, Progress] Merge all annotation information of all the input annotation BED files')
        chroms = [f'chr{n}' for n in range(1, 23)]
        chrom_to_bed_path = {chrom: out_bed_path.replace('.bed', f'.{chrom}.bed') for chrom in chroms}

        # Make a merged BED file for each chromosome
        for chrom in chroms:
            print(f'[{get_curr_time()}, Progress] Merge for {chrom}')
            chrom_merge_bed_path = chrom_to_bed_path[chrom]
            merge_annot_by_chrom(chrom_merge_bed_path, bed_path_dict, chrom)

        # Write headers of the merged bed file
        print(f'[{get_curr_time()}, Progress] Create a BED file with merged annotation information')
        with open(out_bed_path, 'w') as outfile:
            annot_key_str = '|'.join(bed_path_dict.keys())
            print(f'#ANNOT={annot_key_str}', file=outfile)
            print('#chrom', 'start', 'end', 'annot_int', sep='\t', file=outfile)

        # Append the merged BED file of each chromosome
        for chrom in chroms:
            cmd = f'cat {chrom_to_bed_path[chrom]} >> {out_bed_path};'
            cmd += f'rm {chrom_to_bed_path[chrom]};'
            execute_cmd(cmd)

        bgzip_tabix(out_bed_path)


def merge_annot_by_chrom(out_bed_path: str, bed_path_dict: dict, chrom: str):
    """ Merge annotation information of all BED coordinates of one chromosome """
    start_to_key_idx = {}
    end_to_key_idx = {}
    pos_list = []

    # Read coordinates from each annotation bed file
    for i, annot_bed_key in enumerate(bed_path_dict.keys()):
        annot_bed_path = bed_path_dict[annot_bed_key]
        with pysam.TabixFile(annot_bed_path) as bed_file:
            for fields in bed_file.fetch(chrom, parser=pysam.asTuple()):
                start = int(fields[1])
                end = int(fields[2])

                # Init
                if start_to_key_idx.get(start) is None:
                    start_to_key_idx[start] = []

                if end_to_key_idx.get(end) is None:
                    end_to_key_idx[end] = []

                start_to_key_idx[start].append(i)
                end_to_key_idx[end].append(i)
                pos_list.append(start)
                pos_list.append(end)

    # Create a BED file listing new coordinates with merged annotation information
    pos_list.sort()
    one_hot = np.zeros(len(bed_path_dict.keys()))
    prev_pos = -1
    n_bed = 0

    with open(out_bed_path, 'w') as outfile:
        for pos in pos_list:
            if n_bed > 0 and prev_pos != pos:  # Make a new coordinate
                annot_int = one_hot_to_int(one_hot)
                bed_entry = (chrom, prev_pos, pos, annot_int)
                print(*bed_entry, sep='\t', file=outfile)

            end_key_ind = end_to_key_idx.get(pos)
            start_key_ind = start_to_key_idx.get(pos)

            if end_key_ind is not None:  # This position is an end of at least one bed coordinate.
                one_hot[end_key_ind] = 0
                n_bed -= 1

            if start_key_ind is not None:  # This position is a start of at least one bed coordinate.
                one_hot[start_key_ind] = 1
                n_bed += 1

            prev_pos = pos


def one_hot_to_int(one_hot: np.ndarray) -> int:
    n = 0

    for i in range(len(one_hot)):
        if one_hot[i]:
            n += 2 ** i

    return n


if __name__ == '__main__':
    main()

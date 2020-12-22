#!/usr/bin/env python3
"""
A script to concatenate all Z-scores of multiple burden tests
"""
import pandas as pd
from glob import glob
from scipy.stats import binom_test, norm


def main():
    # Results for random mutations
    burden_test_dir = \
        '/data1/mwjeong/projects/cwas/data/' \
        'random-mutation-201010-vep90-cat-burden'
    test_result_paths = sorted(glob(f'{burden_test_dir}/*.txt'))

    # Get a list of original categories
    ori_cat_result_path = \
        '/home/mwjeong/projects/cwas/results/cwas_cat_result_200511.txt'

    with open(ori_cat_result_path) as infile:
        header = infile.readline()
        header_fields = header.strip().split('\t')
        categories = header_fields[1:]

    # Get a default Z-score (Both case and control have values of 0)
    default_p = binom_test(x=1, n=2, p=0.5, alternative='greater')
    default_z = norm.ppf(1 - default_p)

    # Make a DataFrame for Z-scores
    zscore_df = pd.DataFrame()

    for test_result_path in test_result_paths:
        z_dict = get_zscore_dict(test_result_path, default_z, categories)
        zscore_df = zscore_df.append(z_dict, ignore_index=True)

    zscore_df.index.name = 'Simulation'
    zscore_df.index += 1

    # Save the DataFrame as a file
    output_path = \
        '/data1/mwjeong/projects/cwas/data/rand_mut.burden.zscore.201010.txt'
    zscore_df.to_csv(output_path, sep='\t')


def get_zscore_dict(burden_result_path, default_z, categories=None):
    if categories is None:
        z_dict = {}

        with open(burden_result_path) as infile:
            _ = infile.readline()

            for line in infile:
                fields = line.strip().split('\t')
                z_dict[fields[0]] = float(fields[6])
    else:
        z_dict = {category: default_z for category in categories}

        with open(burden_result_path) as infile:
            _ = infile.readline()

            for line in infile:
                fields = line.strip().split('\t')

                if z_dict.get(fields[0]) is not None:
                    z_dict[fields[0]] = float(fields[6])

    return z_dict


if __name__ == '__main__':
    main()

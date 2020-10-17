#!/usr/bin/env python3
"""
A script to estimate and get an effective number of tests
"""
import pickle
import numpy as np
import pandas as pd


# Parse an input z-score matrix
data_dir = '/data1/mwjeong/projects/cwas/data'
zscore_path = f'{data_dir}/rand_mut.burden.zscore.201010.txt'
zscore_df = pd.read_table(zscore_path, index_col='Simulation')
n_col = len(zscore_df.columns.values)
n_row = len(zscore_df.index.values)

# Get a correlation matrix
zscore_mat = zscore_df.values
corr_mat = np.corrcoef(zscore_mat.T)
corr_mat_path = f'{data_dir}/corr_mat_201010.pickle'
pickle.dump(corr_mat, open(corr_mat_path, 'wb'), protocol=4)

# Get a negative laplacian matrix
neg_lap = np.abs(corr_mat)
degrees = np.sum(neg_lap, axis=0)

for i in range(neg_lap.shape[0]):
    neg_lap[i, :] = neg_lap[i, :] / np.sqrt(degrees)
    neg_lap[:, i] = neg_lap[:, i] / np.sqrt(degrees)

neg_lap_path = f'{data_dir}/neg_lap_201010.pickle'
pickle.dump(neg_lap, open(neg_lap_path, 'wb'), protocol=4)

# Get eigen values
eig_vals, _ = np.linalg.eig(neg_lap)
eig_val_path = f'{data_dir}/eig_vals_201010.pickle'
pickle.dump(eig_vals, open(eig_vals, 'wb'))

# Estimate an effective number of tests
e = 1e-12
eig_vals = sorted(eig_vals, key=np.linalg.norm)[::-1]
num_eig_val = min(n_row, n_col)
clean_eig_vals = np.array(eig_vals[:num_eig_val])
clean_eig_vals = clean_eig_vals[clean_eig_vals >= e]
clean_eig_val_total_sum = np.sum(clean_eig_vals)
clean_eig_val_sum = 0
eff_num_test = 0

for i in range(len(clean_eig_vals)):
    clean_eig_val_sum += clean_eig_vals[i]

    if clean_eig_val_sum / clean_eig_val_total_sum >= 0.99:
        eff_num_test = i + 1
        break

print(eff_num_test)

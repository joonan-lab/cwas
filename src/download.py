#!/usr/bin/env python
"""
Script to download essential data for category-wide association study (CWAS)
"""
import os

import yaml

from utils import get_curr_time


def main():
    # Print the script description
    print(__doc__)

    # Paths for this script
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    filepath_conf_path = os.path.join(project_dir, 'conf', 'filepaths.yaml')
    fileurl_conf_path = os.path.join(project_dir, 'conf', 'fileurls.yaml')

    # Parse the configuration files
    with open(filepath_conf_path) as filepath_conf_file:
        filepath_dict = yaml.safe_load(filepath_conf_file)

        for file_group in filepath_dict:
            for file_key in filepath_dict[file_group]:
                filepath_dict[file_group][file_key] = os.path.join(project_dir, filepath_dict[file_group][file_key])

    with open(fileurl_conf_path) as fileurl_conf_file:
        fileurl_dict = yaml.safe_load(fileurl_conf_file)

    # Download the data
    print(f'[{get_curr_time()}, Progress] Download essential data for CWAS')
    download_data(filepath_dict, fileurl_dict)
    print(f'[{get_curr_time()}, Progress] Done')


def download_data(filepath_dict: dict, fileurl_dict: dict):
    """ Download essential data using wget commands"""

    for file_group in fileurl_dict:
        print(f'[{get_curr_time()}, Progress] Start to download data to {file_group}')
        data_dir = filepath_dict[file_group]['data_dir']
        os.makedirs(data_dir, exist_ok=True)

        for data_key in fileurl_dict[file_group]:
            data_dest_path = filepath_dict[file_group][data_key]

            if os.path.isfile(data_dest_path):
                print(f'[{get_curr_time()}, INFO] "{data_dest_path}" already exists. Skip downloading this file.')
            else:
                cmd = f'wget -O {data_dest_path} {fileurl_dict[file_group][data_key]}'
                print(f'[{get_curr_time()}, CMD] {cmd}')
                exit_val = os.system(cmd)

                if exit_val != 0:
                    print(f'[{get_curr_time()}, WARNING] This CMD is failed with this exit value {exit_val}.')


if __name__ == '__main__':
    main()

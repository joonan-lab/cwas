#!/usr/bin/env python
"""
Download essential data for this project
"""
import argparse
import os

import yaml

import cwas.utils.log as log
from cwas.utils.cmd import execute


def main():
    # Print the script description
    print(__doc__)

    # Paths for this script
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(curr_dir)
    filepath_conf_path = \
        os.path.join(project_dir, 'conf', 'download_filepaths.yaml')
    fileurl_conf_path = \
        os.path.join(project_dir, 'conf', 'download_fileurls.yaml')

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--force_overwrite', dest='force_overwrite',
        action='store_const', const=1, default=0,
        help='Force to overwrite when downloading data (Default: 0)'
    )
    args = parser.parse_args()

    # Parse the configuration files
    with open(filepath_conf_path) as filepath_conf_file:
        filepath_dict = yaml.safe_load(filepath_conf_file)

        for file_group in filepath_dict:
            for file_key in filepath_dict[file_group]:
                filepath_dict[file_group][file_key] = \
                    os.path.join(
                        project_dir,
                        filepath_dict[file_group][file_key]
                    )

    with open(fileurl_conf_path) as fileurl_conf_file:
        fileurl_dict = yaml.safe_load(fileurl_conf_file)

    # Download the data
    download_data(filepath_dict, fileurl_dict, args.force_overwrite)
    log.print_progress('Done')


def download_data(filepath_dict: dict, fileurl_dict: dict,
                  force_overwrite: int):
    """ Download essential data using wget commands"""

    for file_group in fileurl_dict:
        log.print_progress(f'Start to download essential data for {file_group}')
        data_dir = filepath_dict[file_group]['data_dir']
        os.makedirs(data_dir, exist_ok=True)

        for data_key in fileurl_dict[file_group]:
            data_dest_path = filepath_dict[file_group][data_key]

            if not force_overwrite and os.path.isfile(data_dest_path):
                log.print_log(
                    'INFO',
                    f'A file "{data_dest_path}" already exists. '
                    f'Skip downloading this file.',
                    True
                )
            else:
                cmd = f'wget -O {data_dest_path} ' \
                      f'{fileurl_dict[file_group][data_key]}'
                execute(cmd)


if __name__ == '__main__':
    main()

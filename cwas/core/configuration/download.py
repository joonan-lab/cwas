"""
Download CWAS default annotation dataset
"""
from cwas.utils.cmd import execute

# Settings
_raw_file_domain = 'https://raw.githubusercontent.com'
_media_file_domain = 'https://media.githubusercontent.com/media'
_username = 'joonan-lab'
_repo_name = 'cwas-dataset'
_branch = 'main'


def download_media_file(dest_dir: str, filename: str):
    remote_filepath = f'{_media_file_domain}/{_username}/' \
                      f'{_repo_name}/{_branch}/{filename}'
    args = ['wget', '-P', dest_dir, remote_filepath]
    execute(args)


def download_raw_file(dest_dir: str, filename: str):
    remote_filepath = f'{_raw_file_domain}/{_username}/' \
                      f'{_repo_name}/{_branch}/{filename}'
    args = ['wget', '-P', dest_dir, remote_filepath]
    execute(args)


def download_bed_tar_gz(dest_dir: str):
    download_media_file(dest_dir, 'bed.tar.gz')


def download_bw_tar_gz(dest_dir: str):
    download_media_file(dest_dir, 'bw.tar.gz')

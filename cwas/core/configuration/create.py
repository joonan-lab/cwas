"""
Create configuration files for CWAS
"""
from __future__ import annotations

import yaml


def create_annotation_key(out_file_path: pathlib.Path,
                          annotation_dir: pathlib.Path, file_type: str):
    assert file_type == 'bed' or file_type == 'bw'

    glob_key = '*.bed.gz' if file_type == 'bed' else '*.bw'
    file_ext_len = len(glob_key) - 1  # Ignore '*'
    filenames = [str(file_path) for file_path in annotation_dir.glob(glob_key)]
    annot_key_dict = {filename: filename[:-file_ext_len].replace('.', '_')
                      for filename in filenames}

    with out_file_path.open('w') as out_f:
        yaml.dump(annot_key_dict, out_f)

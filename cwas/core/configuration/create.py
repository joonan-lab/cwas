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
    filenames = [str(file_path.name)
                 for file_path in annotation_dir.glob(glob_key)]
    annot_key_dict = {filename: filename[:-file_ext_len].replace('.', '_')
                      for filename in filenames}

    with out_file_path.open('w') as out_f:
        yaml.dump(annot_key_dict, out_f)


def split_annotation_key(bed_key_path: pathlib.Path,
                         bw_key_path: pathlib.Path,
                         annotation_key_path: pathlib.Path):
    """Split the input annotation key configuration file"""
    with annotation_key_path.open('r') as annot_key_f:
        annotation_key_dict = yaml.safe_load(annot_key_f)

    bed_key_dict = {}
    bw_key_dict = {}
    for filename in annotation_key_dict:
        if filename.endswith('bed.gz'):
            bed_key_dict[filename] = annotation_key_dict[filename]
        elif filename.endswith('bw'):
            bw_key_dict[filename] = annotation_key_dict[filename]

    with bed_key_path.open('w') as bed_out_f:
        yaml.dump(bed_key_dict, bed_out_f)
    with bw_key_path.open('w') as bw_out_f:
        yaml.dump(bw_key_dict, bw_out_f)

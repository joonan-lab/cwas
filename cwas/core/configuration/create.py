"""
Create configuration files for CWAS
"""
from __future__ import annotations

import yaml

from cwas.core.configuration.settings import get_default_domains


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


def create_category_domain_list(domain_list_path: pathlib.Path,
                                bed_key_path: pathlib.Path,
                                bw_key_path: pathlib.Path,
                                gene_mat_path: pathlib.Path):
    with bed_key_path.open('r') as bed_key_f:
        bed_key_dict = yaml.safe_load(bed_key_f)
        region_domains = bed_key_dict.values()

    with bw_key_path.open('r') as bw_key_f:
        bw_key_dict = yaml.safe_load(bw_key_f)
        cons_domains = bw_key_dict.values()

    with gene_mat_path.open('r') as gene_mat_f:
        header = gene_mat_f.readline()
        columns = header.strip().split()
        gene_list_domains = columns[2:]  # Remove 'gene_id' and 'gene_name'

    domains = get_default_domains()
    domains['region'] += region_domains
    domains['conservation'] += cons_domains
    domains['gene_list'] += gene_list_domains

    with domain_list_path.open('w') as domain_list_f:
        yaml.dump(domains, domain_list_f)

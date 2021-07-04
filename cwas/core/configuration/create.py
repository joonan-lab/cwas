"""
Create configuration files for CWAS
"""
from __future__ import annotations

import yaml

from cwas.core.configuration.settings import get_default_domains
from cwas.core.configuration.settings import get_domain_types
from cwas.core.configuration.settings import get_redundant_domain_pairs
from cwas.utils.log import print_err


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
        yaml.safe_dump(annot_key_dict, out_f)


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
        yaml.safe_dump(bed_key_dict, bed_out_f)
    with bw_key_path.open('w') as bw_out_f:
        yaml.safe_dump(bw_key_dict, bw_out_f)


def create_bw_cutoff_list(bw_cutoff_list: pathlib.Path,
                          bw_key_list: pathlib.Path,
                          user_def_cutoff_list: pathlib.Pathlib = None):
    """ All the arguments must be in a proper YAML format."""
    default_cutoff = 0.0

    with bw_key_list.open('r') as bw_key_f:
        try:
            bw_key_dict = yaml.safe_load(bw_key_f)
            bw_cutoff_dict = {bw_key: default_cutoff
                              for bw_key in bw_key_dict.values()}
        except yaml.YAMLError:
            print_err(f'"{bw_key_list}" is not in a proper YAML format.')
            raise

    with user_def_cutoff_list.open('r') as user_def_cutoff_f:
        try:
            user_def_cutoff_dict = yaml.safe_load(user_def_cutoff_f)
            for bw_filename, bw_cutoff in user_def_cutoff_dict.items():
                bw_key = bw_key_dict.get(bw_filename)
                if bw_key is not None:
                    bw_cutoff_dict[bw_key] = bw_cutoff
        except yaml.YAMLError:
            print_err(f'"{bw_key_list}" is not in a proper YAML format.')
            raise

    with bw_cutoff_list.open('w') as bw_cutoff_f:
        yaml.safe_dump(bw_cutoff_dict, bw_cutoff_f)


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
        yaml.safe_dump(domains, domain_list_f)


def create_redundant_category_table(category_table_path: pathlib.Path):
    domain_types = get_domain_types()
    table_row_template = \
        {domain_type: '*' for domain_type in domain_types}

    with category_table_path.open('w') as out_f:
        print(*domain_types, sep='\t', file=out_f)
        redundant_domain_pairs = get_redundant_domain_pairs()
        for domain_type_pair, domain_pair_set in redundant_domain_pairs.items():
            for domain_pair in domain_pair_set:
                table_row = dict(table_row_template)
                for i, domain_type in enumerate(domain_type_pair):
                    table_row[domain_type] = domain_pair[i]
                print(*[table_row[domain_type] for domain_type in domain_types],
                      sep='\t', file=out_f)

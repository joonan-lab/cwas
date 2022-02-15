"""
Create configuration files for CWAS
"""
from pathlib import Path

import yaml
from cwas.core.configuration.settings import (
    get_default_domains,
    get_domain_types,
    get_redundant_domain_pairs,
)
from cwas.utils.log import print_err


def create_annotation_key(
    out_file_path: Path, annotation_dir: Path, file_type: str
):
    assert file_type == "bed" or file_type == "bw"

    glob_key = "*.bed.gz" if file_type == "bed" else "*.bw"
    file_ext_len = len(glob_key) - 1  # Ignore '*'
    filenames = [
        str(file_path.name) for file_path in annotation_dir.glob(glob_key)
    ]
    annot_key_dict = {
        filename: filename[:-file_ext_len].replace(".", "_")
        for filename in filenames
    }
    _save_as_yaml(out_file_path, annot_key_dict)


def split_annotation_key(
    bed_key_path: Path, bw_key_path: Path, annotation_key_path: Path,
):
    """Split the input annotation key configuration file"""
    annotation_key_dict = _load_yaml_file(annotation_key_path)

    bed_key_dict = {}
    bw_key_dict = {}

    for filename in annotation_key_dict:
        if filename.endswith("bed.gz"):
            bed_key_dict[filename] = annotation_key_dict[filename]
        elif filename.endswith("bw"):
            bw_key_dict[filename] = annotation_key_dict[filename]

    _save_as_yaml(bed_key_path, bed_key_dict)
    _save_as_yaml(bw_key_path, bw_key_dict)


def create_bw_cutoff_list(
    bw_cutoff_list: Path, bw_key_list: Path, user_def_cutoff_list: Path = None,
):
    """ All the arguments must be in a proper YAML format."""
    default_cutoff = 0.0

    bw_key_dict = _load_yaml_file(bw_key_list)
    bw_cutoff_dict = {bw_key: default_cutoff for bw_key in bw_key_dict.values()}
    user_def_cutoff_dict = (
        _load_yaml_file(user_def_cutoff_list)
        if user_def_cutoff_list is not None
        else {}
    )

    for bw_filename, bw_cutoff in user_def_cutoff_dict.items():
        bw_key = bw_key_dict.get(bw_filename)
        if bw_key is not None:
            bw_cutoff_dict[bw_key] = bw_cutoff

    _save_as_yaml(bw_cutoff_list, bw_cutoff_dict)


def create_category_domain_list(
    domain_list_path: Path,
    bed_key_path: Path,
    bw_key_path: Path,
    gene_mat_path: Path,
):
    bed_key_dict = _load_yaml_file(bed_key_path)
    region_domains = bed_key_dict.values()

    bw_key_dict = _load_yaml_file(bw_key_path)
    cons_domains = bw_key_dict.values()

    with gene_mat_path.open("r") as gene_mat_f:
        header = gene_mat_f.readline()
        columns = header.strip().split()
        gene_list_domains = columns[2:]  # Remove 'gene_id' and 'gene_name'

    domains = get_default_domains()
    domains["region"] += region_domains
    domains["conservation"] += cons_domains
    domains["gene_list"] += gene_list_domains

    _save_as_yaml(domain_list_path, domains)


def create_redundant_category_table(category_table_path: Path):
    domain_types = get_domain_types()
    table_row_template = {domain_type: "*" for domain_type in domain_types}

    with category_table_path.open("w") as out_f:
        print(*domain_types, sep="\t", file=out_f)
        redundant_domain_pairs = get_redundant_domain_pairs()
        for domain_type_pair, domain_pair_set in redundant_domain_pairs.items():
            for domain_pair in domain_pair_set:
                table_row = dict(table_row_template)
                for i, domain_type in enumerate(domain_type_pair):
                    table_row[domain_type] = domain_pair[i]
                print(
                    *[table_row[domain_type] for domain_type in domain_types],
                    sep="\t",
                    file=out_f,
                )


def _load_yaml_file(yaml_file_path: Path):
    try:
        with yaml_file_path.open("r") as in_yaml_f:
            yaml_data = yaml.safe_load(in_yaml_f)
    except yaml.YAMLError:
        print_err(f'"{yaml_file_path}" is not in a proper YAML format.')
        raise
    return yaml_data


def _save_as_yaml(yaml_file_path: Path, data):
    with yaml_file_path.open("w") as out_yaml_f:
        yaml.safe_dump(data, out_yaml_f)

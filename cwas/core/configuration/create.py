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

def create_category_domain_list(
    domain_list_path: Path,
    bed_key_path: Path,
    gene_mat_path: Path,
):
    bed_key_dict = _load_yaml_file(bed_key_path)
    
    cons_domains = bed_key_dict['functional_score'].values()
    region_domains = bed_key_dict['functional_annotation'].values()
    
    with gene_mat_path.open("r") as gene_mat_f:
        header = gene_mat_f.readline()
        columns = header.strip().split()
        gene_list_domains = columns[2:]  # Remove 'gene_id' and 'gene_name'

    domains = get_default_domains()
    domains["functional_annotation"] += region_domains
    domains["functional_score"] += cons_domains
    domains["gene_set"] += gene_list_domains

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

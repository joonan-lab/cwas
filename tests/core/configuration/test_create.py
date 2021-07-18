"""
Test cwas.core.configuration.create
"""
import pytest
import yaml

import cwas.core.configuration.create as create


def test_create_annotation_key_bed(cwas_workspace, annotation_dir):
    bed_key_conf = cwas_workspace / 'annotation_key_bed.yaml'
    create.create_annotation_key(bed_key_conf, annotation_dir, 'bed')

    assert bed_key_conf.exists()
    with bed_key_conf.open() as f:
        bed_key = yaml.safe_load(f)
    assert 'bed_annot1.bed.gz' in bed_key
    assert 'bed.annot2.bed.gz' in bed_key
    assert 'bed.annot.bed' not in bed_key
    assert bed_key['bed_annot1.bed.gz'] == 'bed_annot1'
    assert bed_key['bed.annot2.bed.gz'] == 'bed_annot2'

    bed_key_conf.unlink()


def test_create_annotation_key_bw(cwas_workspace, annotation_dir):
    bw_key_conf = cwas_workspace / 'annotation_key_bw.yaml'
    create.create_annotation_key(bw_key_conf, annotation_dir, 'bw')

    assert bw_key_conf.exists()
    with bw_key_conf.open() as f:
        bw_key = yaml.safe_load(f)
    assert 'bw_annot1.bw' in bw_key
    assert 'bw.annot2.bw' in bw_key
    assert 'bw.annot3.bw.gz' not in bw_key
    assert bw_key['bw_annot1.bw'] == 'bw_annot1'
    assert bw_key['bw.annot2.bw'] == 'bw_annot2'

    bw_key_conf.unlink()


def test_split_annotation_key(cwas_workspace, annotation_key_conf):
    bed_key_conf = cwas_workspace / 'split_key_bed.yaml'
    bw_key_conf = cwas_workspace / 'split_key_bw.yaml'
    create.split_annotation_key(bed_key_conf, bw_key_conf, annotation_key_conf)

    assert bed_key_conf.exists()
    with bed_key_conf.open() as f:
        bed_key = yaml.safe_load(f)
    assert 'bed_annot1.bed.gz' in bed_key and bed_key['bed_annot1.bed.gz'] == \
           'bed1'
    assert 'bed.annot2.bed.gz' in bed_key and bed_key['bed.annot2.bed.gz'] == \
           'bed2'
    assert 'bed.annot3.bed' not in bed_key

    assert bw_key_conf.exists()
    with bw_key_conf.open() as f:
        bw_key = yaml.safe_load(f)
    assert 'bw_annot1.bw' in bw_key and bw_key['bw_annot1.bw'] == 'bw1'
    assert 'bw.annot2.bw' in bw_key and bw_key['bw.annot2.bw'] == 'bw2'
    assert 'bw.annot3.bw.gz' not in bw_key

    bed_key_conf.unlink()
    bw_key_conf.unlink()


def test_create_bw_cutoff_list(cwas_workspace, annotation_dir, bw_cutoff_conf):
    bw_cutoff_list = cwas_workspace / 'annotation_cutoff_bw.yaml'
    bw_key_list = cwas_workspace / 'annotation_key_bw.yaml'
    create.create_annotation_key(bw_key_list, annotation_dir, 'bw')
    create.create_bw_cutoff_list(bw_cutoff_list, bw_key_list, bw_cutoff_conf)

    assert bw_cutoff_list.exists()

    with bw_key_list.open() as f:
        bw_key_dict = yaml.safe_load(f)
    with bw_cutoff_list.open() as f:
        bw_cutoff_dict = yaml.safe_load(f)
    with bw_cutoff_conf.open() as f:
        user_def_cutoff_dict = yaml.safe_load(f)

    for bw_filename, bw_key in bw_key_dict.items():
        assert bw_key in bw_cutoff_dict and \
               user_def_cutoff_dict[bw_filename] == bw_cutoff_dict[bw_key]
    assert not (set(bw_cutoff_dict.keys()) - set(bw_key_dict.values()))

    bw_cutoff_list.unlink()
    bw_key_list.unlink()


def test_create_bw_cutoff_wo_user_def(cwas_workspace, annotation_dir):
    bw_cutoff_list = cwas_workspace / 'annotation_cutoff_bw.yaml'
    bw_key_list = cwas_workspace / 'annotation_key_bw.yaml'
    create.create_annotation_key(bw_key_list, annotation_dir, 'bw')
    create.create_bw_cutoff_list(bw_cutoff_list, bw_key_list)

    assert bw_cutoff_list.exists()

    with bw_key_list.open() as f:
        bw_key_dict = yaml.safe_load(f)
    with bw_cutoff_list.open() as f:
        bw_cutoff_dict = yaml.safe_load(f)
    assert not (set(bw_cutoff_dict.keys()) - set(bw_key_dict.values()))

    bw_cutoff_list.unlink()
    bw_key_list.unlink()


def test_create_category_domain_list(cwas_workspace, annotation_key_conf,
                                     gene_matrix):
    # Setting
    bed_key_conf = cwas_workspace / 'split_key_bed.yaml'
    bw_key_conf = cwas_workspace / 'split_key_bw.yaml'
    create.split_annotation_key(bed_key_conf, bw_key_conf, annotation_key_conf)

    domain_list_path = cwas_workspace / 'category_domain.yaml'
    create.create_category_domain_list(domain_list_path, bed_key_conf,
                                       bw_key_conf, gene_matrix)

    assert domain_list_path.is_file()

    with domain_list_path.open('r') as domain_list_f, \
            bed_key_conf.open('r') as bed_key_f, \
            bw_key_conf.open('r') as bw_key_f, \
            gene_matrix.open('r') as gene_mat_f:
        domain_dict = yaml.safe_load(domain_list_f)
        bed_key_dict = yaml.safe_load(bed_key_f)
        bw_key_dict = yaml.safe_load(bw_key_f)
        gene_mat_header = gene_mat_f.readline()
        gene_list_domains = gene_mat_header.strip().split()
        gene_list_domains = gene_list_domains[2:]

    assert len(domain_dict['region']) == 1 + len(bed_key_dict)
    assert len(domain_dict['conservation']) == 1 + len(bw_key_dict)
    assert len(domain_dict['gene_list']) == 1 + len(gene_list_domains)

    bed_key_conf.unlink()
    bw_key_conf.unlink()
    domain_list_path.unlink()


def test_create_redundant_category_table(cwas_workspace):
    redundant_category_table = cwas_workspace / 'redundant_category.txt'
    create.create_redundant_category_table(redundant_category_table)
    assert redundant_category_table.exists()
    redundant_category_table.unlink()


def test__load_yaml_file_error(cwas_workspace):
    test_path = cwas_workspace / 'not_yaml.txt'
    with test_path.open('w') as test_f:
        print('{Hello: World!}"', file=test_f)
    with pytest.raises(yaml.YAMLError):
        _ = create.split_annotation_key(None, None, test_path)
    test_path.unlink()

"""
Test cwas.configuration
"""
from cwas.env import Env
from cwas.configuration import Configuration


def test_run_configuration(cwas_workspace, annotation_dir, annotation_key_conf,
                           bw_cutoff_conf, gene_matrix):
    argv = [
        '-d', str(annotation_dir),
        '-m', str(gene_matrix),
        '-k', str(annotation_key_conf),
        '-c', str(bw_cutoff_conf),
        '-w', str(cwas_workspace),
    ]
    inst = Configuration.get_instance(argv)
    setattr(inst, 'env', Env())
    inst.run()

    data_dir = cwas_workspace / 'annotation-data'
    gene_matrix = cwas_workspace / 'gene_matrix.txt'
    bed_key_list = cwas_workspace / 'annotation_key_bed.yaml'
    bw_key_list = cwas_workspace / 'annotation_key_bw.yaml'
    bw_cutoff_list = cwas_workspace / 'annotation_cutoff_bw.yaml'
    category_domain_list = cwas_workspace / 'category_domain.yaml'
    redundant_category_table = cwas_workspace / 'redundant_category.txt'

    assert data_dir.is_dir()
    assert gene_matrix.is_file()
    assert bed_key_list.is_file()
    assert bw_key_list.is_file()
    assert bw_cutoff_list.is_file()
    assert category_domain_list.is_file()
    assert redundant_category_table.is_file()

    # Teardown
    for f in data_dir.glob('*'):
        f.unlink()
    data_dir.unlink()
    gene_matrix.unlink()
    bed_key_list.unlink()
    bw_key_list.unlink()
    bw_cutoff_list.unlink()
    category_domain_list.unlink()
    redundant_category_table.unlink()

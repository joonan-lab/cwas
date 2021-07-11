"""
Test cwas.configuration
"""
from cwas.env import Env
from cwas.configuration import Configuration


def test_run_configuration(tmp_dir, annotation_dir, annotation_key_conf,
                           bw_cutoff_conf, gene_matrix):
    argv = [
        '-d', str(annotation_dir),
        '-m', str(gene_matrix),
        '-k', str(annotation_key_conf),
        '-c', str(bw_cutoff_conf),
        '-w', str(tmp_dir),
    ]
    inst = Configuration.get_instance(argv)
    setattr(inst, 'env', Env())
    inst.run()

    data_dir = tmp_dir / 'annotation-data'
    gene_matrix = tmp_dir / 'gene_matrix.txt'
    bed_key_list = tmp_dir / 'annotation_key_bed.yaml'
    bw_key_list = tmp_dir / 'annotation_key_bw.yaml'
    bw_cutoff_list = tmp_dir / 'annotation_cutoff_bw.yaml'
    category_domain_list = tmp_dir / 'category_domain.yaml'
    redundant_category_table = tmp_dir / 'redundant_category.txt'

    assert data_dir.is_dir()
    assert gene_matrix.is_file()
    assert bed_key_list.is_file()
    assert bw_key_list.is_file()
    assert bw_cutoff_list.is_file()
    assert category_domain_list.is_file()
    assert redundant_category_table.is_file()

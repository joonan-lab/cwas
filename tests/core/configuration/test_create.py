"""
Test cwas.core.configuration.create
"""
import cwas.core.configuration.create as create
import pytest
import yaml


@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace):
    cwas_workspace.mkdir()


@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace):
    yield
    remove_workspace(cwas_workspace)


def remove_workspace(cwas_workspace):
    for f in cwas_workspace.glob("*"):
        f.unlink()
    cwas_workspace.rmdir()


def test_create_annotation_key_bed(cwas_workspace, annotation_dir):
    bed_key_conf = cwas_workspace / "annotation_keys.yaml"
    create.create_annotation_key(bed_key_conf, annotation_dir, "bed")

    assert bed_key_conf.exists()
    with bed_key_conf.open() as f:
        bed_key = yaml.safe_load(f)
    assert "bed_annot1.bed.gz" in bed_key
    assert "bed.annot2.bed.gz" in bed_key
    assert "bed.annot.bed" not in bed_key
    assert bed_key["bed_annot1.bed.gz"] == "bed_annot1"
    assert bed_key["bed.annot2.bed.gz"] == "bed_annot2"

    bed_key_conf.unlink()


def test_create_category_domain_list(
    cwas_workspace, annotation_key_conf, gene_matrix
):
    # Setting
    bed_key_conf = cwas_workspace / "annotation_keys.yaml"
    create.split_annotation_key(bed_key_conf, annotation_key_conf)

    domain_list_path = cwas_workspace / "category_domain.yaml"
    create.create_category_domain_list(
        domain_list_path, bed_key_conf, gene_matrix
    )

    assert domain_list_path.is_file()

    with domain_list_path.open("r") as domain_list_f, bed_key_conf.open(
        "r"
    ) as bed_key_f, gene_matrix.open(
        "r"
    ) as gene_mat_f:
        domain_dict = yaml.safe_load(domain_list_f)
        bed_key_dict = yaml.safe_load(bed_key_f)
        region_dict = bed_key_dict['functional_annotation']
        score_dict = bed_key_dict['functional_score']
        gene_mat_header = gene_mat_f.readline()
        gene_list_domains = gene_mat_header.strip().split()
        gene_list_domains = gene_list_domains[2:]

    assert len(domain_dict["region"]) == 1 + len(region_dict)
    assert len(domain_dict["conservation"]) == 1 + len(score_dict)
    assert len(domain_dict["gene_list"]) == 1 + len(gene_list_domains)

    bed_key_conf.unlink()
    domain_list_path.unlink()


def test_create_redundant_category_table(cwas_workspace):
    redundant_category_table = cwas_workspace / "redundant_category.txt"
    create.create_redundant_category_table(redundant_category_table)
    assert redundant_category_table.exists()
    redundant_category_table.unlink()


def test__load_yaml_file_error(cwas_workspace):
    test_path = cwas_workspace / "not_yaml.txt"
    with test_path.open("w") as test_f:
        print('{Hello: World!}"', file=test_f)
    with pytest.raises(yaml.YAMLError):
        _ = create.split_annotation_key(None, None, test_path)
    test_path.unlink()

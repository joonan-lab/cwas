"""
Test cwas.preparation
"""
import pytest
import yaml
from cwas.annotation import Annotation
from cwas.env import Env
import cwas.cli
import sys


class AnnotationMock(Annotation):
    """Mocking the Annotation class"""

    def run(self):
        pass


@pytest.fixture(scope="module")
def vcf_path(cwas_workspace):
    vcf = cwas_workspace / "test_target.vcf"
    create_vcf_file(vcf)
    return vcf

@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace, annotation_dir):
    cwas_workspace.mkdir(exist_ok=True)  # This line creates the workspace directory
    set_env(cwas_workspace, annotation_dir)
    annotation_key_path = cwas_workspace / "annotation_keys.yaml"
    annotation_key_path.parent.mkdir(parents=True, exist_ok=True)  # Ensure parent directory exists
    create_annotation_key_yaml(annotation_key_path)

@pytest.fixture(scope="module")
def output_dir_path(cwas_workspace):
    output_dir = cwas_workspace 
    return output_dir

@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace):
    yield
    env = Env()
    env.reset()
    env.remove_file()
    for f in cwas_workspace.glob("*"):
        f.unlink()
    cwas_workspace.rmdir()


def set_env(cwas_workspace, annotation_dir):
    env = Env()
    env.set_env("CWAS_WORKSPACE", cwas_workspace)
    env.set_env("VEP", "VEP")
    env.set_env("ANNOTATION_DATA", annotation_dir)
    env.set_env("ANNOTATION_BED_KEY", cwas_workspace / "annotation_keys.yaml")
    env.set_env("MERGED_BED", cwas_workspace / "merged.bed.gz")
    env.set_env("MERGED_BED_INDEX", cwas_workspace / "merged.bed.gz.tbi")
    env.save()

def create_annotation_key_yaml(annot_key_path):
    key_data = {
        "functional_score": {
            "bed_annot1.bed.gz": "bed1",
            "bed.annot2.bed": "bed2"
        },
        "functional_annotation": {
            "bed.annot3.bed": "bed3",
            "bed_annot4.bed.gz": "bed4"
        }
    }
    with annot_key_path.open("w") as outfile:
        yaml.safe_dump(key_data, outfile) 

def create_vcf_file(vcf_path):
    vcf_header = (
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
    )
    vcf_entries = [
        ("chr1", "100", ".", "A", "G", ".", "PASS", "."),
        ("chr1", "200", ".", "C", "T", ".", "PASS", "."),
        ("chr2", "300", ".", "G", "A", ".", "PASS", "."),
    ]
    with vcf_path.open("w") as vcf_file:
        print(*vcf_header, sep="\t", file=vcf_file)
        for entry in vcf_entries:
            print(*entry, sep="\t", file=vcf_file)

@pytest.fixture(scope="module")
def required_args(vcf_path, output_dir_path):
    return ["-v", str(vcf_path), "-o", str(output_dir_path)]


def test_parse_args(required_args, vcf_path, output_dir_path):
    sys.argv = ['cwas', 'annotation', *required_args]
    inst = cwas.cli.main()
    assert getattr(inst.args, "vcf_path") == vcf_path
    assert getattr(inst.args, "output_dir_path") == output_dir_path

def test_parse_args_without_required_arg():
    with pytest.raises(SystemExit):
        sys.argv = ['cwas', 'annotation']
        cwas.cli.main()


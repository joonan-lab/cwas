"""
Test cwas.preparation
"""
import random
from multiprocessing import cpu_count

import pytest
import yaml
from cwas.annotation import Annotation
from cwas.env import Env


class AnnotationMock(Annotation):
    """Mocking the Annotation class"""

    def run(self):
        pass


@pytest.fixture(scope="module")
def vcf_path(cwas_workspace):
    return cwas_workspace / "test_target.vcf"


@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace, annotation_dir, vcf_path):
    cwas_workspace.mkdir()
    set_env(cwas_workspace, annotation_dir)
    create_vcf_file(vcf_path)
    create_bw_key_yaml(annotation_dir / "bw_key.yaml")


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
    env.set_env("ANNOTATION_BW_KEY", annotation_dir / "bw_key.yaml")
    env.set_env("MERGED_BED", cwas_workspace / "merged.bed.gz")
    env.set_env("MERGED_BED_INDEX", cwas_workspace / "merged.bed.gz.tbi")
    env.save()


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
    # TODO: Add VCF entries (variants)
    with vcf_path.open("w") as vcf_file:
        print(*vcf_header, sep="\t", file=vcf_file)


def create_bw_key_yaml(bw_key_yaml_path):
    bw_key_data = {f"test{i + 1}.bw": f"test{i + 1}" for i in range(3)}
    with bw_key_yaml_path.open("w") as outfile:
        yaml.safe_dump(bw_key_data, outfile)


@pytest.fixture(scope="module")
def required_args(vcf_path):
    return ["-v", str(vcf_path)]


def test_parse_args(required_args, vcf_path):
    cpu = random.choice(range(1, cpu_count() + 1))
    args = required_args + ["-p", str(cpu)]
    inst = AnnotationMock.get_instance(args)
    assert getattr(inst, "vcf_path") == vcf_path
    assert getattr(inst, "num_proc") == cpu


def test_parse_args_required_arg_only(required_args, vcf_path):
    inst = AnnotationMock.get_instance(required_args)
    assert getattr(inst, "vcf_path") == vcf_path
    assert getattr(inst, "num_proc") == 1


def test_parse_args_without_required_arg():
    args = []
    with pytest.raises(SystemExit):
        AnnotationMock.get_instance(args)


def test_parse_args_invalid_num_proc(required_args):
    cpu = cpu_count() + 1
    args = required_args + ["-p", str(cpu)]
    with pytest.raises(ValueError):
        AnnotationMock.get_instance(args)

    args = required_args + ["-p", "0"]
    with pytest.raises(ValueError):
        AnnotationMock.get_instance(args)

    args = required_args + ["-p", "-1"]
    with pytest.raises(ValueError):
        AnnotationMock.get_instance(args)


def test_get_bw_custom_annotations(required_args, annotation_dir):
    inst = AnnotationMock.get_instance(required_args)
    for i, bw_custom_annotation in enumerate(inst.bw_custom_annotations):
        assert (
            str(annotation_dir / f"test{i + 1}.bw"),
            f"test{i + 1}",
        ) == bw_custom_annotation

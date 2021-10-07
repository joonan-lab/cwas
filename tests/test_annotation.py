"""
Test cwas.preparation
"""
import random
from multiprocessing import cpu_count

import pytest
from cwas.env import Env
from cwas.annotation import Annotation


class AnnotationWithoutRun(Annotation):
    """Mocking the Annotation class"""

    def run(self):
        pass


@pytest.fixture(scope="module")
def env_mock(cwas_workspace):
    env = Env()
    env.set_env("CWAS_WORKSPACE", cwas_workspace)
    env.set_env("VEP", "VEP")
    return env


def test_run_annotation_with_empty_env():
    inst = Annotation.get_instance([])
    setattr(inst, "env", Env())
    with pytest.raises(RuntimeError):
        inst.run()


def test_parse_args():
    vcf_path = "test.vcf"
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-i", vcf_path, "-p", str(cpu), "--per_chrom"]
    inst = AnnotationWithoutRun.get_instance(args)
    assert getattr(inst, "vcf_path") == vcf_path
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "per_chrom") == 1


def test_parse_args_without_num_proc():
    vcf_path = "test.vcf"
    args = ["-i", vcf_path, "--per_chrom"]
    inst = AnnotationWithoutRun.get_instance(args)
    assert getattr(inst, "vcf_path") == vcf_path
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "per_chrom") == 1


def test_parse_args_required_only():
    vcf_path = "test.vcf"
    args = ["-i", vcf_path]
    inst = AnnotationWithoutRun.get_instance(args)
    assert getattr(inst, "vcf_path") == vcf_path
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "per_chrom") == 0


def test_parse_args_without_required_arg():
    args = []
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)


def test_parse_args_invalid_num_proc():
    vcf_path = "test.vcf"
    base_args = ["-i", vcf_path, "--per_chrom"]

    cpu = cpu_count() + 1
    args = base_args + ["-p", str(cpu)]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

    args = base_args + ["-p", "0"]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

    args = base_args + ["-p", "-1"]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)


def test_parse_args_num_proc_without_per_chrom():
    # '-p [int]' argument cannot be used without '--per_chrom' argument
    vcf_path = "test.vcf"
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-i", vcf_path, "-p", str(cpu)]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

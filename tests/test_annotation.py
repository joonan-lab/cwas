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


def test_default_args():
    inst = AnnotationWithoutRun.get_instance([])
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "per_chrom") == 0


def test_parse_args():
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-p", str(cpu), "--per_chrom"]
    inst = AnnotationWithoutRun.get_instance(args)
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "per_chrom") == 1

    args = ["--per_chrom"]
    inst = AnnotationWithoutRun.get_instance(args)
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "per_chrom") == 1


def test_parse_args_value_error():
    # '-p [int]' argument cannot be used without '--per_chrom' argument
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-p", str(cpu)]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)
    
    cpu = cpu_count() + 1
    args = ["-p", str(cpu)]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

    args = ["-p", "0"]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

    args = ["-p", "-1"]
    with pytest.raises(ValueError):
        AnnotationWithoutRun.get_instance(args)

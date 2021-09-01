"""
Test cwas.preparation
"""
import random
from multiprocessing import cpu_count

import pytest
from cwas.env import Env
from cwas.preparation import Preparation


class PreparationWithoutRun(Preparation):
    """Mocking the Preparation class"""

    def run(self):
        pass


@pytest.fixture(scope="module")
def env_mock(cwas_workspace, annotation_dir):
    env = Env()
    env.set_env("CWAS_WORKSPACE", cwas_workspace)
    env.set_env("ANNOTATION_DATA", annotation_dir)
    env.set_env("ANNOTATION_BED_KEY", cwas_workspace / "bed_key.yaml")
    return env


def test_run_preparation_with_empty_env():
    inst = Preparation.get_instance([])
    setattr(inst, "env", Env())
    with pytest.raises(RuntimeError):
        inst.run()


def test_default_args():
    inst = PreparationWithoutRun.get_instance([])
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "force_overwrite") == 0


def test_parse_args():
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-p", str(cpu), "--force_overwrite"]
    inst = PreparationWithoutRun.get_instance(args)
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "force_overwrite") == 1

    args = ["-p", str(cpu)]
    inst = PreparationWithoutRun.get_instance(args)
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "force_overwrite") == 0

    args = ["--force_overwrite"]
    inst = PreparationWithoutRun.get_instance(args)
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "force_overwrite") == 1


def test_parse_args_value_error():
    cpu = cpu_count() + 1
    args = ["-p", str(cpu)]
    with pytest.raises(ValueError):
        PreparationWithoutRun.get_instance(args)

    args = ["-p", "0"]
    with pytest.raises(ValueError):
        PreparationWithoutRun.get_instance(args)

    args = ["-p", "-1"]
    with pytest.raises(ValueError):
        PreparationWithoutRun.get_instance(args)

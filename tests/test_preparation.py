"""
Test cwas.preparation
"""
import random
from multiprocessing import cpu_count
from pathlib import Path
from typing import Tuple

import pytest
from cwas.env import Env
from cwas.preparation import Preparation


class PreparationMock(Preparation):
    """Mocking the Preparation class"""

    def _prepare_annotation(self) -> Tuple[Path, Path]:
        """This step is supposed to be tested in 
        tests/core/preparation/test_annotation.py
        """
        merged_bed_path = self.workspace / "merged.bed.gz"
        merged_bed_idx_path = self.workspace / "merged.bed.gz.tbi"
        return merged_bed_path, merged_bed_idx_path


@pytest.fixture(scope="function")
def set_env(cwas_workspace, annotation_dir):
    env = Env()
    env_path = cwas_workspace / ".cwas_config"
    env.set_path(env_path)
    env.set_env("CWAS_WORKSPACE", cwas_workspace)
    env.set_env("ANNOTATION_DATA", annotation_dir)
    env.set_env("ANNOTATION_BED_KEY", cwas_workspace / "bed_key.yaml")
    env.save()


def test_run_without_configuration():
    inst = PreparationMock.get_instance()
    assert inst.get_env("CWAS_WORKSPACE") is None
    assert inst.get_env("ANNOTATION_DATA") is None
    assert inst.get_env("ANNOTATION_BED_KEY") is None
    with pytest.raises(RuntimeError):
        inst.run()


def test_default_args():
    inst = PreparationMock.get_instance()
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "force_overwrite") == 0


def test_parse_args():
    cpu = random.choice(range(1, cpu_count() + 1))
    args = ["-p", str(cpu), "--force_overwrite"]
    inst = PreparationMock.get_instance(args)
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "force_overwrite") == 1

    args = ["-p", str(cpu)]
    inst = PreparationMock.get_instance(args)
    assert getattr(inst, "num_proc") == cpu
    assert getattr(inst, "force_overwrite") == 0

    args = ["--force_overwrite"]
    inst = PreparationMock.get_instance(args)
    assert getattr(inst, "num_proc") == 1
    assert getattr(inst, "force_overwrite") == 1


def test_parse_args_value_error():
    cpu = cpu_count() + 1
    args = ["-p", str(cpu)]
    with pytest.raises(ValueError):
        PreparationMock.get_instance(args)

    args = ["-p", "0"]
    with pytest.raises(ValueError):
        PreparationMock.get_instance(args)

    args = ["-p", "-1"]
    with pytest.raises(ValueError):
        PreparationMock.get_instance(args)


def test_env_after_run_preparation(set_env):
    inst = PreparationMock.get_instance()
    inst.run()
    assert inst.get_env("MERGED_BED")
    assert inst.get_env("MERGED_BED_INDEX")

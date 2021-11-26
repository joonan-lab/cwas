import shutil

import pytest
from cwas.core.annotation.vep import VEP


@pytest.fixture(scope="module")
def installed_vep() -> str:
    return shutil.which("vep")


def test_init_vep(installed_vep):
    vep_inst = VEP(installed_vep)
    assert vep_inst.get_vep_path() == installed_vep


def test_init_vep_with_invalid_path(cwas_workspace):
    invalid_vep_path = str(cwas_workspace / "_vep_not_exists")
    with pytest.raises(ValueError):
        _ = VEP(invalid_vep_path)


def test_init_vep_with_no_arg():
    with pytest.raises(ValueError):
        _ = VEP()


def test_init_vep_with_none():
    with pytest.raises(ValueError):
        _ = VEP(None)

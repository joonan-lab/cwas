import shutil
from pathlib import Path

import pytest
from cwas.core.annotation.vep import VEP


@pytest.fixture(scope="module")
def installed_vep() -> str:
    return shutil.which("vep")


@pytest.fixture(scope="module")
def input_vcf_path(cwas_workspace) -> Path:
    tmp_vcf_path = cwas_workspace / "test_vep.vcf"
    tmp_vcf_path.touch()
    yield tmp_vcf_path
    tmp_vcf_path.unlink()


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


def test_set_input_vcf(installed_vep, input_vcf_path):
    vep_inst = VEP(installed_vep)
    vep_inst.set_input_vcf(input_vcf_path)
    assert f"-i {input_vcf_path}" in vep_inst.get_cmd()


def test_set_invalid_input_vcf(installed_vep, cwas_workspace):
    vep_inst = VEP(installed_vep)
    invalid_vcf_path = cwas_workspace / "test_not_exists.vcf"
    with pytest.raises(FileNotFoundError):
        vep_inst.set_input_vcf(invalid_vcf_path)


def test_set_input_vcf_with_none():
    vep_inst = VEP(installed_vep)
    with pytest.raises(ValueError):
        vep_inst.set_input_vcf(None)

import shutil
from pathlib import Path

import pytest
from cwas.core.annotation.vep import VEP


@pytest.fixture(scope="module")
def installed_vep_path() -> str:
    return shutil.which("vep")


@pytest.fixture(scope="module")
def input_vcf_path(cwas_workspace) -> Path:
    return cwas_workspace / "test_vep.vcf"


@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace, input_vcf_path):
    cwas_workspace.mkdir()
    input_vcf_path.touch()


@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace, input_vcf_path):
    yield
    input_vcf_path.unlink()
    cwas_workspace.rmdir()


def test_init_vep(installed_vep_path):
    vep_inst = VEP(installed_vep_path)
    assert vep_inst.get_vep_path() == installed_vep_path


def test_init_vep_with_invalid_path(cwas_workspace):
    invalid_vep_path = str(cwas_workspace / "_vep_not_exists")
    with pytest.raises(FileNotFoundError):
        _ = VEP(invalid_vep_path)


def test_init_vep_with_no_arg():
    with pytest.raises(TypeError):
        _ = VEP()


def test_init_vep_with_none():
    with pytest.raises(ValueError):
        _ = VEP(None)


def test_set_input_vcf(installed_vep_path, input_vcf_path):
    vep_inst = VEP(installed_vep_path)
    vep_inst.set_input_vcf(str(input_vcf_path))
    assert f"-i {input_vcf_path}" in vep_inst.get_cmd()


def test_set_invalid_input_vcf(installed_vep_path, cwas_workspace):
    vep_inst = VEP(installed_vep_path)
    invalid_vcf_path = str(cwas_workspace / "test_not_exists.vcf")
    with pytest.raises(FileNotFoundError):
        vep_inst.set_input_vcf(invalid_vcf_path)


def test_set_input_vcf_with_none(installed_vep_path):
    vep_inst = VEP(installed_vep_path)
    with pytest.raises(ValueError):
        vep_inst.set_input_vcf(None)


def test_output_vcf_path(installed_vep_path, input_vcf_path):
    vep_inst = VEP(installed_vep_path)
    vep_inst.set_input_vcf(str(input_vcf_path))
    assert vep_inst.output_vcf_path == str(input_vcf_path).replace(
        ".vcf", ".annotated.vcf"
    )


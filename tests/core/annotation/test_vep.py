from pathlib import Path

import pytest
from cwas.core.annotation.vep import VepCmdGenerator


@pytest.fixture(scope="module")
def vep_path(cwas_workspace) -> str:
    return str(cwas_workspace / "vep_mock")

@pytest.fixture(scope="module")
def input_vcf_path(cwas_workspace) -> str:
    return str(cwas_workspace / "test_vep.vcf")

@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace, vep_path, input_vcf_path):
    cwas_workspace.mkdir()
    Path(vep_path).touch()
    Path(input_vcf_path).touch()

@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace, vep_path, input_vcf_path):
    yield
    Path(vep_path).unlink()
    Path(input_vcf_path).unlink()
    cwas_workspace.rmdir()


def test_init_vep(vep_path, input_vcf_path,
                  vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    vep_inst = VepCmdGenerator(vep_path=vep_path, input_vcf_path=input_vcf_path,
                               vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                               vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                               vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)
    assert vep_inst.vep_path == vep_path
    assert vep_inst.input_vcf_path == input_vcf_path


def test_init_vep_with_invalid_vep_path(cwas_workspace, input_vcf_path, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    invalid_vep_path = str(cwas_workspace / "_vep_not_exists")
    with pytest.raises(FileNotFoundError):
        _ = VepCmdGenerator(vep_path=invalid_vep_path, input_vcf_path=input_vcf_path,
                            vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                            vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                            vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)


def test_init_vep_with_no_arg():
    with pytest.raises(TypeError):
        _ = VepCmdGenerator()


def test_init_vep_with_no_vep(input_vcf_path, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    with pytest.raises(ValueError):
        _ = VepCmdGenerator(vep_path=None, input_vcf_path=input_vcf_path,
                            vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                            vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                            vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)


def test_init_vep_with_no_input_vcf(vep_path, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    with pytest.raises(ValueError):
        _ = VepCmdGenerator(vep_path=vep_path, input_vcf_path=None,
                            vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                            vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                            vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)


def test_init_vep_with_invalid_input_vcf(vep_path, cwas_workspace, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    invalid_vcf_path = str(cwas_workspace / "test_not_exists.vcf")
    with pytest.raises(FileNotFoundError):
        _ = VepCmdGenerator(vep_path=vep_path, input_vcf_path=invalid_vcf_path,
                            vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                            vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                            vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)


def test_output_vcf_path(vep_path, input_vcf_path, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    vep_inst = VepCmdGenerator(vep_path=vep_path, input_vcf_path=input_vcf_path,
                               vep_cache_path=vep_dir, vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                               vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                               vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)
    assert vep_inst.output_vcf_path == input_vcf_path.replace(
        ".vcf", ".vep.vcf"
    )


def test_cmd(vep_path, input_vcf_path, vep_dir, vep_conserv, vep_loftee, vep_ances, vep_gerp, vep_msdb, vep_mskey):
    vep_inst = VepCmdGenerator(vep_path=vep_path, input_vcf_path=input_vcf_path,
                               vep_cache_path=str(vep_dir), vep_conservation_path=vep_conserv, vep_loftee_path=vep_loftee,
                               vep_human_ancestor_fa_path=vep_ances, vep_gerp_bw_path=vep_gerp,
                               vep_mis_db_path=vep_msdb, vep_mis_info_key=vep_mskey, num_proc=1)
    assert vep_inst.cmd_str.startswith(vep_path)
    assert f"-i {input_vcf_path}" in vep_inst.cmd_str
    assert (
        f"-o {input_vcf_path.replace('.vcf', '.vep.vcf')}" in vep_inst.cmd_str
    )

@pytest.fixture(scope="module", autouse=True)
def create_vep_dir(
    vep_dir, vep_conserv, vep_loftee, vep_msdb, vep_ances, vep_gerp
):
    vep_dir.mkdir()
    #create_vep_conserv(vep_conserv)
    #create_misdb(vep_msdb)
    Path(vep_conserv).touch()
    Path(vep_msdb).touch()
    Path(vep_loftee).mkdir()
    Path(vep_ances).touch()
    Path(vep_gerp).touch()
    print("[TEST] Temporary VEP directory has created.")
    yield
    Path(vep_loftee).rmdir()
    for f in vep_dir.glob("*"):
        f.unlink()
    vep_dir.rmdir()
    print("[TEST] Temporary VEP directory has deleted.")

#def test_cmd_for_bw_custom_annotation(vep_path, input_vcf_path, annotation_dir):
#    vep_inst = VepCmdGenerator(vep_path, input_vcf_path)
#    bw_paths = [
#        (bw_path, f"test{i + 1}")
#        for i, bw_path in enumerate(annotation_dir.glob("*.bw"))
#    ]
#    for bw_path, bw_key in bw_paths:
#        vep_inst.add_bw_custom_annotation(str(bw_path), bw_key)

#    for bw_path, bw_key in bw_paths:
#        assert (
#            f"--custom {bw_path},{bw_key},bigwig,overlap,0" in vep_inst.cmd_str
#        )

from pathlib import Path

from cwas.utils.check import check_is_file


class VEP:
    def __init__(self, vep: str) -> None:
        self.vep = vep
        self._check_vep_path()
        self.input_vcf_path = None

    def get_vep_path(self):
        return self.vep

    def set_input_vcf(self, input_vcf_path: Path):
        self.input_vcf_path = input_vcf_path
        self._check_input_vcf_path()

    def _check_vep_path(self):
        check_is_file(self.vep)

    def _check_input_vcf_path(self):
        check_is_file(self.input_vcf_path)

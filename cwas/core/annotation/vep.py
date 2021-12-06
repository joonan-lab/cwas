from cwas.utils.check import check_is_file


class VEP:
    def __init__(self, vep_path: str) -> None:
        self.vep_path = vep_path
        self._check_vep_path()
        self.input_vcf_path = None

    def get_vep_path(self) -> str:
        return self.vep_path

    def set_input_vcf(self, input_vcf_path: str):
        self.input_vcf_path = input_vcf_path
        self._check_input_vcf_path()

    def get_cmd(self) -> str:
        return " ".join([self.vep_path, "-i", str(self.input_vcf_path)])

    def _check_vep_path(self):
        check_is_file(self.vep_path)

    def _check_input_vcf_path(self):
        check_is_file(self.input_vcf_path)

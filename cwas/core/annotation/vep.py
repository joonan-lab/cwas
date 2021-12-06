from cwas.utils.check import check_is_file


class VEP:
    def __init__(self, vep_path: str, input_vcf_path: str) -> None:
        check_is_file(vep_path)
        check_is_file(input_vcf_path)
        self._vep_path = vep_path
        self._input_vcf_path = input_vcf_path
        self._output_vcf_path = input_vcf_path.replace(".vcf", ".annotated.vcf")

    def get_vep_path(self) -> str:
        return self._vep_path

    def get_input_vcf_path(self) -> str:
        return self._input_vcf_path

    def get_output_vcf_path(self) -> str:
        return self._output_vcf_path

    def get_cmd(self) -> str:
        return " ".join(
            [
                self._vep_path,
                "-i",
                self._input_vcf_path,
                "-o",
                self._output_vcf_path,
            ]
        )

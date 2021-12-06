from cwas.utils.check import check_is_file


class VEP:
    def __init__(self, vep_path: str, input_vcf_path: str) -> None:
        check_is_file(vep_path)
        check_is_file(input_vcf_path)
        self.vep_path = vep_path
        self.input_vcf_path = input_vcf_path
        self.output_vcf_path = input_vcf_path.replace(".vcf", ".annotated.vcf")

    def get_vep_path(self) -> str:
        return self.vep_path

    def get_cmd(self) -> str:
        return " ".join([self.vep_path, "-i", str(self.input_vcf_path)])

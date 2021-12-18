"""
Command line generator for Variant Effect Predictor (VEP)
"""
from cwas.utils.check import check_is_file


class VepCmdGenerator:
    def __init__(self, vep_path: str, input_vcf_path: str) -> None:
        self._vep_path = vep_path
        self._check_vep_path()
        self._input_vcf_path = input_vcf_path
        self._check_input_vcf_path()
        self._output_vcf_path = input_vcf_path.replace(".vcf", ".annotated.vcf")

    def _check_vep_path(self):
        try:
            check_is_file(self._vep_path)
        except ValueError:
            raise ValueError(f"Invalid VEP path: {self._vep_path}")
        except:
            raise

    def _check_input_vcf_path(self):
        try:
            check_is_file(self._input_vcf_path)
        except ValueError:
            raise ValueError(f"Invalid VCF path: {self._input_vcf_path}")
        except:
            raise

    @property
    def vep_path(self) -> str:
        return self._vep_path

    @property
    def input_vcf_path(self) -> str:
        return self._input_vcf_path

    @property
    def output_vcf_path(self) -> str:
        return self._output_vcf_path

    @property
    def cmd(self) -> str:
        args = [
            self._vep_path,
            "-i",
            self._input_vcf_path,
            "-o",
            self._output_vcf_path,
        ]
        args += self._get_cmd_basic_vep_opt()
        args += self._get_cmd_opt_pick_one_gene_isoform()
        args += self._get_cmd_opt_pick_nearest_gene()
        return " ".join(args)

    def _get_cmd_basic_vep_opt(self) -> list:
        """Return basic options (no plugins) of VEP"""
        return [
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--format",
            "vcf",
            "--vcf",
            "--no_stats",
            "--polyphen p",
        ]

    def _get_cmd_opt_pick_one_gene_isoform(self) -> list:
        """Return options in order to pick a gene isoform 
        with most severe consequence"""
        return [
            "--per_gene",
            "--pick",
            "--pick_order",
            "canonical,appris,tsl,biotype,ccds,rank,length",
        ]

    def _get_cmd_opt_pick_nearest_gene(self) -> list:
        """Return options in order to pick the nearest gene"""
        return ["--distance", "2000", "--nearest", "symbol", "--symbol"]

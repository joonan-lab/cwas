"""
Command line generator for Variant Effect Predictor (VEP)
"""
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir


class VepCmdGenerator:
    def __init__(self, vep_path: str,
                 vep_conservation_path: str, vep_loftee_path: str, vep_human_ancestor_fa_path: str, vep_gerp_bw_path: str, vep_mpc_path: str,
                 input_vcf_path: str) -> None:
        self._vep_path = vep_path
        self._vep_conservation_path = vep_conservation_path
        self._vep_loftee_path = vep_loftee_path
        self._vep_human_ancestor_fa_path = vep_human_ancestor_fa_path
        self._vep_gerp_bw_path = vep_gerp_bw_path
        self._vep_mpc_path = vep_mpc_path
        self._check_vep_path()
        self._check_vep_conservation_path()
        self._check_vep_loftee_path()
        self._check_vep_human_ancestor_fa_path()
        self._check_vep_gerp_bw_path()
        self._check_vep_mpc_path()
        self._input_vcf_path = input_vcf_path
        self._check_input_vcf_path()
        self._output_vcf_path = input_vcf_path.replace(".vcf", ".vep.vcf")
        self._bw_custom_annotations = []

    def _check_vep_path(self):
        try:
            check_is_file(self._vep_path)
        except ValueError:
            raise ValueError(f"Invalid VEP path: {self._vep_path}")
        except Exception:
            raise

    def _check_vep_conservation_path(self):
        try:
            check_is_file(self._vep_conservation_path)
        except ValueError:
            raise ValueError(f"Invalid VEP resource path (conservation file): {self._vep_conservation_path}")
        except Exception:
            raise

    def _check_vep_loftee_path(self):
        try:
            check_is_dir(self._vep_loftee_path)
        except ValueError:
            raise ValueError(f"Invalid VEP resource path (loftee directory): {self._vep_loftee_path}")
        except Exception:
            raise

    def _check_vep_human_ancestor_fa_path(self):
        try:
            check_is_file(self._vep_human_ancestor_fa_path)
        except ValueError:
            raise ValueError(f"Invalid VEP resource path (human ancestor): {self._vep_human_ancestor_fa_path}")
        except Exception:
            raise

    def _check_vep_gerp_bw_path(self):
        try:
            check_is_file(self._vep_gerp_bw_path)
        except ValueError:
            raise ValueError(f"Invalid VEP resource path (gerp bigwig): {self._vep_gerp_bw_path}")
        except Exception:
            raise

    def _check_vep_mpc_path(self):
        try:
            check_is_file(self._vep_mpc_path)
        except ValueError:
            raise ValueError(f"Invalid VEP resource path (MPC): {self._vep_mpc_path}")
        except Exception:
            raise

    def _check_input_vcf_path(self):
        try:
            check_is_file(self._input_vcf_path)
        except ValueError:
            raise ValueError(f"Invalid VCF path: {self._input_vcf_path}")
        except Exception:
            raise

    def add_bw_custom_annotation(self, bw_path: str, annotation_key: str):
        check_is_file(bw_path)
        self._bw_custom_annotations.append((bw_path, annotation_key))

    @property
    def vep_path(self) -> str:
        return self._vep_path

    @property
    def vep_conservation_path(self) -> str:
        return self._vep_conservation_path

    @property
    def vep_loftee_path(self) -> str:
        return self._vep_loftee_path

    @property
    def vep_human_ancestor_fa_path(self) -> str:
        return self._vep_human_ancestor_fa_path

    @property
    def vep_gerp_bw_path(self) -> str:
        return self._vep_gerp_bw_path

    @property
    def vep_mpc_path(self) -> str:
        return self._vep_mpc_path

    @property
    def input_vcf_path(self) -> str:
        return self._input_vcf_path

    @property
    def output_vcf_path(self) -> str:
        return self._output_vcf_path

    @output_vcf_path.setter
    def output_vcf_path(self, arg: str):
        self._output_vcf_path = arg

    @property
    def cmd_str(self) -> str:
        return " ".join(self.cmd)

    @property
    def cmd(self) -> list:
        result = [
            self._vep_path,
            "-i",
            self._input_vcf_path,
            "-o",
            self._output_vcf_path,
        ]
        result += self.cmd_option_basic
        result += self.cmd_option_pick_one_gene_isoform
        result += self.cmd_option_pick_nearest_gene
        result += self.cmd_option_bw_custom_annotations
        return result

    @property
    def cmd_option_basic(self) -> list:
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
            "--plugin",
            ''.join(['LoF,conservation_file:', self._vep_conservation_path,
                     ',loftee_path:', self._vep_loftee_path,
                     ',human_ancestor_fa:', self._vep_human_ancestor_fa_path,
                     '/human_ancestor.fa.gz,gerp_bigwig:', self._vep_gerp_bw_path]),
            "--dir_plugins",
            self._vep_loftee_path,
            "--plugin",
            ''.join(["MPC,", self._vep_mpc_path]),
        ]

    @property
    def cmd_option_pick_one_gene_isoform(self) -> list:
        """Return options in order to pick a gene isoform
        with most severe consequence"""
        return [
            "--per_gene",
            "--pick",
            "--pick_order",
            "rank,canonical,appris,tsl,biotype,ccds,length",
        ]

    @property
    def cmd_option_pick_nearest_gene(self) -> list:
        """Return options in order to pick the nearest gene"""
        return ["--distance", "2000", "--nearest", "symbol", "--symbol"]

    @property
    def cmd_option_bw_custom_annotations(self) -> list:
        result = []

        for bw_path, annotation_key in self._bw_custom_annotations:
            result += [
                "--custom",
                ",".join([bw_path, annotation_key, "bigwig", "overlap", "0"]),
            ]

        return result

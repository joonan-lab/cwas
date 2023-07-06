"""
Command line generator for Variant Effect Predictor (VEP)
"""
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir


class VepCmdGenerator:
    def __init__(self, vep_path: str,
                 vep_cache_path: str, vep_conservation_path: str,
                 vep_loftee_path: str, vep_human_ancestor_fa_path: str,
                 vep_gerp_bw_path: str, vep_mis_db_path: str, 
                 vep_mis_info_key: str, input_vcf_path: str, 
                 num_proc: str) -> None:
        self._vep_path = vep_path
        self._vep_cache_path = vep_cache_path
        self._vep_conservation_path = vep_conservation_path
        self._vep_loftee_path = vep_loftee_path
        self._vep_human_ancestor_fa_path = vep_human_ancestor_fa_path
        self._vep_gerp_bw_path = vep_gerp_bw_path
        self._vep_mis_db_path = vep_mis_db_path
        self._vep_mis_info_key = vep_mis_info_key
        self._input_vcf_path = input_vcf_path
        self._check_validity()
        self._output_vcf_path = input_vcf_path.replace(".vcf", ".vep.vcf")
        self._num_proc = num_proc

    @staticmethod
    def _check_path(path: str, message: str, is_dir: bool = False):
        try:
             if is_dir:
                 check_is_dir(path)
             else:
                 check_is_file(path)
        except ValueError:
            raise ValueError(f"{message}: {path}")
        except Exception:
            raise

    def _check_validity(self):
         self._check_path(self._vep_path, "Invalid VEP path")
         self._check_path(self._vep_conservation_path, "Invalid VEP resource path (conservation file)")
         self._check_path(self._vep_loftee_path, "Invalid VEP resource path (loftee directory)", is_dir=True)
         self._check_path(self._vep_human_ancestor_fa_path, "Invalid VEP resource path (human ancestor fasta file)")
         self._check_path(self._vep_gerp_bw_path, "Invalid VEP resource path (gerp bigwig file)")
         self._check_path(self._vep_mis_db_path, "Invalid VEP resource path (missense database file)")
         self._check_path(self._input_vcf_path, "Invalid input VCF path")
         self._check_path(self._vep_cache_path, "Invalid VEP cache directory path", is_dir=True)

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
    def vep_mis_db_path(self) -> str:
        return self._vep_mis_db_path

    @property
    def vep_mis_info_key(self) -> str:
        return self._vep_mis_info_key

    @property
    def vep_cache_path(self) -> str:
        return self._vep_cache_path

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
    def num_proc(self) -> str:
        return self._num_proc

    @property
    def cmd_str(self) -> str:
        return " ".join(self.cmd)

    @property
    def cmd(self) -> list:
        result = [
            self.vep_path,
            "-i",
            self.input_vcf_path,
            "-o",
            self.output_vcf_path,
        ]
        result += self.cmd_option_basic
        result += self.cmd_option_pick_one_gene_isoform
        result += self.cmd_option_pick_nearest_gene
        return result

    @property
    def cmd_option_basic(self) -> list:
        """Return basic options (no plugins) of VEP"""
        return [
            "--assembly",
            "GRCh38",
            "--offline",
            "--cache",
            "--dir_cache",
            self.vep_cache_path,
            "--force_overwrite",
            "--format",
            "vcf",
            "--vcf",
            "--no_stats",
            "--plugin",
            ','.join(['LoF',
                      'conservation_file:' + self.vep_conservation_path,
                      'loftee_path:' + self.vep_loftee_path,
                      'human_ancestor_fa:' + self.vep_human_ancestor_fa_path,
                      'gerp_bigwig:' + self.vep_gerp_bw_path]),
            "--dir_plugins",
            self.vep_loftee_path,
            "--custom",
            ",".join([self._vep_mis_db_path, "MisDb", 'vcf', "exact", "0", self.vep_mis_info_key]),
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
import argparse
import shutil
from pathlib import Path

import cwas.core.configuration.create as create
import cwas.utils.log as log
import cwas.utils.check as check
from cwas.runnable import Runnable


class Configuration(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self.workspace = self._get_workspace()
        self.user_config = self.workspace / "configuration.txt"
        self.data_dir_symlink = self.workspace / "annotation-data"
        self.gene_matrix_symlink = self.workspace / "gene_matrix.txt"
        self.bed_key_list = self.workspace / "annotation_key_bed.yaml"
        self.bw_key_list = self.workspace / "annotation_key_bw.yaml"
        self.bw_cutoff_list = self.workspace / "annotation_cutoff_bw.yaml"
        self.category_domain_list = self.workspace / "category_domain.yaml"
        self.redundant_category_table = (
            self.workspace / "redundant_category.txt"
        )

    def _get_workspace(self):
        workspace = self.get_env("CWAS_WORKSPACE")
        if workspace is None:
            raise RuntimeError(
                "${CWAS_WORKSPACE} does not set. Run the Start step first."
            )
        return Path(workspace)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        """This argparse is only supposed to print help
        in order to explain what each configuration means.
        """
        parser = argparse.ArgumentParser(
            description="Arguments for CWAS Configuration",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument(
            "-d",
            "--annotation_data_dir",
            dest="data_dir",
            required=False,
            type=Path,
            help="Path to your annotation data directory",
        )
        parser.add_argument(
            "-m",
            "--gene_matrix",
            dest="gene_matrix",
            required=False,
            type=Path,
            help="Path to your gene matrix",
        )
        parser.add_argument(
            "-a",
            "--annotation_key_config",
            dest="annot_key_conf",
            required=False,
            type=Path,
            help="Path to a configuration file (.yaml) that "
            "specifies the annotation key of each "
            "annotation data file",
        )
        parser.add_argument(
            "-b",
            "--bigwig_cutoff_config",
            dest="bw_cutoff_conf",
            required=False,
            type=Path,
            help="Path to a configuration file (.yaml) that "
            "specifies the annotation cutoff of "
            "each BigWig file",
        )
        parser.add_argument(
            "-v",
            "--vep",
            dest="vep",
            required=False,
            type=Path,
            help="Path to Variant Effect Predictor (VEP)",
        )
        parser.add_argument(
            "-vconv",
            "--vep_conservation",
            dest="vep_conservation",
            required=False,
            type=Path,
            help="Path to your VEP resource (conservation file)",
        )
        parser.add_argument(
            "-vlof",
            "--vep_loftee",
            dest="vep_loftee",
            required=False,
            type=Path,
            help="Path to your VEP resource (loftee directory)",
        )
        parser.add_argument(
            "-vha",
            "--vep_human_ancestor_fa",
            dest="vep_human_ancestor_fa",
            required=False,
            type=Path,
            help="Path to your VEP resource (human ancestor)",
        )
        parser.add_argument(
            "-vgerp",
            "--vep_gerp_bw",
            dest="vep_gerp_bw",
            required=False,
            type=Path,
            help="Path to your VEP resource (gerp bigwig)",
        )
        parser.add_argument(
            "-vmpc",
            "--vep_mpc",
            dest="vep_mpc",
            required=False,
            type=Path,
            help="Path to your VEP resource (MPC)",
        )
        return parser

    def run(self):
        self._set_config_to_attr()
        self._check_attr_from_user_config()
        self._create_data_dir_symlink()
        self._create_gene_matrix_symlink()
        self._create_annotation_key_list()
        self._create_bw_cutoff_list()
        self._create_category_info()
        self._set_env()

    def _set_config_to_attr(self):
        user_config = self._load_configuration()
        self.data_dir = Path(user_config.get("ANNOTATION_DATA_DIR"))
        self.gene_matrix = Path(user_config.get("GENE_MATRIX"))
        self.vep = Path(user_config.get("VEP"))
        self.vep_conservation = Path(user_config.get("VEP_CONSERVATION_FILE"))
        self.vep_loftee = Path(user_config.get("VEP_LOFTEE"))
        self.vep_human_ancestor_fa = Path(user_config.get("VEP_HUMAN_ANCESTOR_FA"))
        self.vep_gerp_bw = Path(user_config.get("VEP_GERP_BIGWIG"))
        self.vep_mpc = Path(user_config.get("VEP_MPC"))

        annot_key_conf = user_config.get("ANNOTATION_KEY_CONFIG")
        bw_cutoff_conf = user_config.get("BIGWIG_CUTOFF_CONFIG")
        self.annot_key_conf = (
            None if not annot_key_conf else Path(annot_key_conf)
        )
        self.bw_cutoff_conf = (
            None if not bw_cutoff_conf else Path(bw_cutoff_conf)
        )

    def _load_configuration(self) -> dict:
        user_config_dict = {}

        with self.user_config.open() as user_config_file:
            for line in user_config_file:
                k, v = line.strip().split("=")
                # The two 'CONFIG' file paths are optional.
                if not k.endswith("CONFIG"):
                    self._check_config_value(k, v)
                user_config_dict[k] = v

        return user_config_dict

    def _check_config_value(self, key: str, value: str):
        if not value:
            raise ValueError(f"The value for '{key}' is empty.")

    def _check_attr_from_user_config(self):
        check.check_is_file(self.user_config)
        check.check_is_dir(self.data_dir)
        check.check_is_file(self.vep_conservation)
        check.check_is_dir(self.vep_loftee)
        check.check_is_file(self.vep_human_ancestor_fa)
        check.check_is_file(self.vep_gerp_bw)
        check.check_is_file(self.vep_mpc)
        check.check_is_file(self.gene_matrix)
        if self.annot_key_conf is not None:
            check.check_is_file(self.annot_key_conf)
        if self.bw_cutoff_conf is not None:
            check.check_is_file(self.bw_cutoff_conf)

        if self.vep is not None and shutil.which(self.vep) is None:
            raise ValueError(f'"{self.vep} is not an executable."')

    def _create_data_dir_symlink(self):
        data_dir = getattr(self, "data_dir")
        data_dir_symlink = getattr(self, "data_dir_symlink")
        log.print_progress(
            f"Create a symlink to your annotation data "
            f'directory "{data_dir_symlink}"'
        )
        try:
            data_dir_symlink.symlink_to(data_dir, target_is_directory=True)
        except FileExistsError:
            log.print_warn(
                f'"{data_dir_symlink}" already exists so skip '
                f"creating the symbolic link"
            )

    def _create_gene_matrix_symlink(self):
        gene_matrix = getattr(self, "gene_matrix")
        gene_matrix_symlink = getattr(self, "gene_matrix_symlink")
        log.print_progress(
            f"Create a symlink to your gene matrix " f'"{gene_matrix_symlink}"'
        )

        try:
            gene_matrix_symlink.symlink_to(gene_matrix)
        except FileExistsError:
            log.print_warn(
                f'"{gene_matrix_symlink}" already exists so skip '
                f"creating the symbolic link"
            )

    def _create_annotation_key_list(self):
        bed_key_list = getattr(self, "bed_key_list")
        bw_key_list = getattr(self, "bw_key_list")
        log.print_progress(
            f"Create annotation key lists "
            f'"{bed_key_list}" and "{bw_key_list}"'
        )
        data_dir = getattr(self, "data_dir_symlink")
        annot_key_conf = getattr(self, "annot_key_conf")

        if annot_key_conf is None:
            create.create_annotation_key(bed_key_list, data_dir, "bed")
            create.create_annotation_key(bw_key_list, data_dir, "bw")
        else:
            create.split_annotation_key(
                bed_key_list, bw_key_list, annot_key_conf
            )

    def _create_bw_cutoff_list(self):
        bw_cutoff_list = getattr(self, "bw_cutoff_list")
        bw_key_list = getattr(self, "bw_key_list")
        bw_cutoff_conf = getattr(self, "bw_cutoff_conf")  # User-defined
        log.print_progress(f'Create BigWig cufoff list "{bw_cutoff_list}"')
        create.create_bw_cutoff_list(
            bw_cutoff_list, bw_key_list, bw_cutoff_conf
        )

    def _create_category_info(self):
        """ Create a list of category domains and a redundant category table"""
        domain_list = getattr(self, "category_domain_list")
        bed_key_list = getattr(self, "bed_key_list")
        bw_key_list = getattr(self, "bw_key_list")
        gene_matrix = getattr(self, "gene_matrix")
        redundant_category_table = getattr(self, "redundant_category_table")

        log.print_progress(
            f"Create a CWAS category domain list " f'"{domain_list}"'
        )
        create.create_category_domain_list(
            domain_list, bed_key_list, bw_key_list, gene_matrix
        )
        log.print_progress(
            f"Create a redundant category table "
            f'"{redundant_category_table}"'
        )
        create.create_redundant_category_table(redundant_category_table)

    def _set_env(self):
        log.print_progress("Set CWAS environment variables")
        self.set_env("VEP", getattr(self, "vep"))
        self.set_env("VEP_CONSERVATION_FILE", getattr(self, "vep_conservation"))
        self.set_env("VEP_LOFTEE", getattr(self, "vep_loftee"))
        self.set_env("VEP_HUMAN_ANCESTOR_FA", getattr(self, "vep_human_ancestor_fa"))
        self.set_env("VEP_GERP_BIGWIG", getattr(self, "vep_gerp_bw"))
        self.set_env("VEP_MPC", getattr(self, "vep_mpc"))
        self.set_env("ANNOTATION_DATA", self.data_dir_symlink)
        self.set_env("GENE_MATRIX", self.gene_matrix_symlink)
        self.set_env("ANNOTATION_BED_KEY", self.bed_key_list)
        self.set_env("ANNOTATION_BW_KEY", self.bw_key_list)
        self.set_env("ANNOTATION_BW_CUTOFF", self.bw_cutoff_list)
        self.set_env("CATEGORY_DOMAIN", self.category_domain_list)
        self.set_env("REDUNDANT_CATEGORY", self.redundant_category_table)
        self.save_env()

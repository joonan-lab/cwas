"""
CWAS Annotation Step

This step annotate user's VCF file using annotation data specified
in the CWAS configuration step. This step mainly uses
Variant Effect Predictor (VEP) to annotate user's VCF file.
"""
import argparse
from pathlib import Path

import yaml

from cwas.core.annotation.bed import annotate as _annotate_using_bed
from cwas.core.annotation.vep import VepCmdGenerator
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir
from cwas.utils.cmd import CmdExecutor, compress_using_bgzip, index_using_tabix
from cwas.utils.log import print_arg, print_log, print_progress


class Annotation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments of CWAS annotation step",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        default_workspace = Path.home() / ".cwas"
        parser.add_argument(
            "-v",
            "--vcf_file",
            dest="vcf_path",
            required=True,
            type=Path,
            help="Target VCF file",
        )
        parser.add_argument(
            "-o_dir",
            "--output_directory",
            dest="output_dir_path",
            required=False,
            default=default_workspace,
            type=Path,
            help="Directory where output file will be saved",
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Target VCF file", args.vcf_path)
        print_arg("Target VCF file", args.output_dir_path)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_file(args.vcf_path)
        check_is_dir(args.output_dir_path)

    @property
    def vcf_path(self):
        return self.args.vcf_path.resolve()

    @property
    def output_dir_path(self):
        return self.args.output_dir_path.resolve()

    @property
    def vep_cmd(self):
        vep_cmd_generator = VepCmdGenerator(
            self.get_env("VEP"),
            self.get_env("VEP_CONSERVATION_FILE"), self.get_env("VEP_LOFTEE"), self.get_env("VEP_HUMAN_ANCESTOR_FA"), self.get_env("VEP_GERP_BIGWIG"), self.get_env("VEP_MPC"),
            str(self.vcf_path)
        )
        vep_cmd_generator.output_vcf_path = self.vep_output_vcf_path
        for bw_path, annotation_key in self.bw_custom_annotations:
            vep_cmd_generator.add_bw_custom_annotation(bw_path, annotation_key)
        return vep_cmd_generator.cmd

    @property
    def vep_output_vcf_path(self):
        return (
            f"{self.output_dir_path}/"
            f"{self.vcf_path.name.replace('.vcf', '.vep.vcf')}"
        )

    @property
    def vep_output_vcf_gz_path(self):
        return self.vep_output_vcf_path.replace(".vcf", ".vcf.gz")

    @property
    def annotated_vcf_path(self):
        return (
            f"{self.output_dir_path}/"
            f"{self.vcf_path.name.replace('.vcf', '.annotated.vcf')}"
        )

    @property
    def bw_custom_annotations(self):
        with open(self.get_env("ANNOTATION_BW_KEY")) as infile:
            bw_custom_path_dict = yaml.safe_load(infile)
        annotation_data_dir = self.get_env("ANNOTATION_DATA")

        for bw_filename, bw_annotation_key in bw_custom_path_dict.items():
            yield (f"{annotation_data_dir}/{bw_filename}", bw_annotation_key)

    def run(self):
        self.annotate_using_bigwig()
        self.process_vep_vcf()
        self.annotate_using_bed()
        self.update_env()

    def annotate_using_bigwig(self):
        print_progress("BigWig custom annotations via VEP")
        if (
            Path(self.vep_output_vcf_path).is_file()
            or Path(self.vep_output_vcf_gz_path).is_file()
        ):
            print_log(
                "NOTICE",
                "You have already done the BigWig custom annotations.",
                True,
            )
            return

        vep_bin, *vep_args = self.vep_cmd
        CmdExecutor(vep_bin, vep_args).execute_raising_err()

    def process_vep_vcf(self):
        print_progress("Compress the VEP output using bgzip")
        vcf_gz_path = compress_using_bgzip(self.vep_output_vcf_path)
        print_progress("Create an index of the VEP output using tabix")
        index_using_tabix(vcf_gz_path)

    def annotate_using_bed(self):
        print_progress("BED custom annotation")
        if Path(self.annotated_vcf_path).is_file():
            print_log(
                "NOTICE",
                "You have already done the BED custom annotation.",
                True,
            )
            return
        _annotate_using_bed(
            self.vep_output_vcf_gz_path,
            self.annotated_vcf_path,
            self.get_env("MERGED_BED"),
        )

    def update_env(self):
        self.set_env("ANNOTATED_VCF", self.annotated_vcf_path)
        self.save_env()

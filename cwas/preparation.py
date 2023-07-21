import argparse
from pathlib import Path
from typing import Tuple

import yaml

import cwas.utils.log as log
from cwas.core.preparation.annotation import merge_bed_files
from cwas.runnable import Runnable
from cwas.utils.check import check_num_proc
from cwas.utils.cmd import compress_using_bgzip, index_using_tabix


class Preparation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg("No. Processes for this step", args.num_proc)
        log.print_arg(
            "Force to overwrite the result",
            "Y" if args.force_overwrite else "N",
        )

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)

    @property
    def num_proc(self):
        return self.args.num_proc

    @property
    def force_overwrite(self):
        return bool(self.args.force_overwrite)

    def run(self):
        self._load_env()
        bed_gz_path, bed_idx_path = self._prepare_annotation()
        self._save_as_env(bed_gz_path, bed_idx_path)

    def _load_env(self):
        """Load environment variables to attributes"""
        try:
            self.workspace = Path(self.get_env("CWAS_WORKSPACE"))
            self.annot_data_dir = Path(self.get_env("ANNOTATION_DATA"))
            self.bed_key_list_path = self.workspace / self.get_env(
                "ANNOTATION_BED_KEY"
            )
        except TypeError:
            raise RuntimeError(
                "Failed to get one of CWAS environment variable."
                " Maybe you omitted to run Configuration step."
            )

    def _prepare_annotation(self) -> Tuple[Path, Path]:
        log.print_progress(
            "Data preprocessing to prepare CWAS annotation step"
        )

        with self.bed_key_list_path.open() as bed_key_list_file:
            bed_key_list = yaml.safe_load(bed_key_list_file)
        
        # Merge two dictionaries
        bed_key_list = bed_key_list['functional_score'] | bed_key_list['functional_annotation']

        bed_file_and_keys = []
        for bed_filename, bed_key in bed_key_list.items():
            bed_file_path = self.annot_data_dir / bed_filename
            bed_file_and_keys.append((bed_file_path, bed_key))

        log.print_progress(
            "Merge all of your annotation BED files into one BED file"
        )
        merge_bed_path = self.workspace / "merged_annotation.bed"
        merge_bed_files(
            merge_bed_path,
            bed_file_and_keys,
            self.num_proc,
            self.force_overwrite,
        )
        log.print_progress("Compress your BED file.")
        bed_gz_path = compress_using_bgzip(
            merge_bed_path, self.force_overwrite
        )

        log.print_progress("Make an index of your BED file.")
        bed_idx_path = index_using_tabix(bed_gz_path, self.force_overwrite)

        return bed_gz_path, bed_idx_path

    def _save_as_env(self, bed_gz_path: Path, bed_idx_path: Path):
        self.set_env("MERGED_BED", bed_gz_path)
        self.set_env("MERGED_BED_INDEX", bed_idx_path)
        self.save_env()

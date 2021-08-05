import argparse
from pathlib import Path

import yaml

import cwas.utils.log as log
from cwas.core.preparation.annotation import (compress_bed_file,
                                              index_bed_file, merge_bed_files)
from cwas.runnable import Runnable
from cwas.utils.error import check_num_proc


class Preparation(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description='Arguments for Annotation Data Preparation',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument('-p', '--num_proc', dest='num_proc',
                            required=False, type=int, default=1,
                            help='Max No. processes for this step')
        parser.add_argument('-f', '--force_overwrite', dest='force_overwrite',
                            action='store_const', const=1, default=0,
                            help='Force to overwrite the result')

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg('No. Processes for this step', args.num_proc)
        log.print_arg('Force to overwrite the result',
                      'Y' if args.force_overwrite else 'N')

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)

    def run(self):
        self._prepare_annotation()

    def _prepare_annotation(self):
        log.print_progress(
            'Data preprocessing to prepare CWAS annotation step')
        cwas_env = getattr(self, 'env')
        workspace = Path(cwas_env.get_env('CWAS_WORKSPACE'))
        annot_data_dir = Path(cwas_env.get_env('ANNOTATION_DATA'))
        bed_key_list_path = workspace / cwas_env.get_env('ANNOTATION_BED_KEY')

        with bed_key_list_path.open() as bed_key_list_file:
            bed_key_list = yaml.safe_load(bed_key_list_file)

        bed_file_and_keys = []
        for bed_filename, bed_key in bed_key_list.items():
            bed_file_path = annot_data_dir / bed_filename
            bed_file_and_keys.append((bed_file_path, bed_key))

        log.print_progress(
            'Merge all of your annotation BED files into one BED file')
        num_proc = getattr(self, 'num_proc')
        force_overwrite = getattr(self, 'force_overwrite')
        merge_bed_path = workspace / 'merged_annotation.bed'
        merge_bed_files(merge_bed_path, bed_file_and_keys,
                        num_proc, force_overwrite)

        log.print_progress("Compress your BED file.")
        bed_gz_path = compress_bed_file(merge_bed_path)

        log.print_progress("Make an index of your BED file.")
        bed_idx_path = index_bed_file(bed_gz_path)

        cwas_env.set_env('MERGED_BED', bed_gz_path)
        cwas_env.set_env('MERGED_BED_INDEX', bed_idx_path)

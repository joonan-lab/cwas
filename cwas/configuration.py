import argparse
from pathlib import Path

import cwas.core.configuration.create as create
import cwas.utils.log as log
from cwas.env import set_env
from cwas.runnable import Runnable


class Configuration(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        work_dir = getattr(self, 'work_dir')
        self.data_dir_symlink = work_dir / 'annotation-data'
        self.gene_mat_symlink = work_dir / 'gene_matrix.txt'
        self.bed_key_conf = work_dir / 'annotation_key_bed.yaml'
        self.bw_key_conf = work_dir / 'annotation_key_bw.yaml'

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description='Arguments for CWAS Configuration',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        default_work_dir = Path.home() / '.cwas'
        parser.add_argument('-d', '--annotation_data_dir', dest='data_dir',
                            required=True, type=Path,
                            help="Path to your annotation data directory")
        parser.add_argument('-m', '--gene_matrix', dest='gene_matrix',
                            required=True, type=Path,
                            help="Path to your gene matrix")
        parser.add_argument('-k', '--annotation_key', dest='annot_key_conf',
                            required=False, type=Path, default=None,
                            help="Path to a configuration file (.yaml) that "
                                 "specifies the annotation key of each "
                                 "annotation data file")
        parser.add_argument('-c', '--bigwig_cutoff', dest='bw_cutoff_conf',
                            required=False, type=Path, default=None,
                            help="Path to a configuration file (.yaml) that "
                                 "specifies the annotation cutoff of "
                                 "each BigWig file")
        parser.add_argument('-w', '--workspace', dest='work_dir',
                            required=False,
                            type=Path, default=default_work_dir,
                            help='Path to your CWAS workspace')
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg('Your annotation data directory', args.data_dir)
        log.print_arg('Your gene matrix', args.gene_matrix)
        log.print_arg('Your annotation key list', args.annot_key_conf)
        log.print_arg('Your BigWig cutoff list', args.bw_cutoff_conf)
        log.print_arg('CWAS workspace', args.work_dir)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        if args.work_dir.exists():
            if args.work_dir.is_dir():
                log.print_warn('The CWAS workspace already exists.')
            else:
                log.print_err('The argument points to another kind of file.')
                raise NotADirectoryError(
                    f"Non-directory file: '{args.work_dir}'")

    def run(self):
        self._create_workspace()
        set_env('CWAS_WORKSPACE', str(getattr(self, 'work_dir')))
        self._create_data_dir_symlink()
        self._create_gene_matrix_symlink()
        self._create_annotation_key_list()

        log.print_log('Notice', 'Not implemented yet.')

    def _create_workspace(self):
        work_dir = getattr(self, 'work_dir')
        log.print_progress(f'Create CWAS workspace "{work_dir}"')
        try:
            work_dir.mkdir(parents=True, exist_ok=True)
        except NotADirectoryError:
            log.print_err('The path to CWAS workspace is invalid.')
            raise

    def _create_data_dir_symlink(self):
        data_dir = getattr(self, 'data_dir')
        data_dir_symlink = getattr(self, 'data_dir_symlink')
        log.print_progress(f'Create a symlink to your annotation data '
                           f'directory "{data_dir_symlink}"')
        try:
            data_dir_symlink.symlink_to(data_dir, target_is_directory=True)
        except FileExistsError:
            log.print_warn(f'"{data_dir_symlink}" already exists so skip '
                           f'creating the symbolic link')

    def _create_gene_matrix_symlink(self):
        gene_matrix = getattr(self, 'gene_matrix')
        gene_matrix_symlink = getattr(self, 'gene_matrix_symlink')
        log.print_progress(f'Create a symlink to your gene matrix '
                           f'"{gene_matrix_symlink}"')

        try:
            gene_matrix_symlink.symlink_to(gene_matrix)
        except FileExistsError:
            log.print_warn(f'"{gene_matrix_symlink}" already exists so skip '
                           f'creating the symbolic link')


    def _create_annotation_key_list(self):
        bed_key_conf = getattr(self, 'bed_key_conf')
        bw_key_conf = getattr(self, 'bw_key_conf')
        log.print_progress(f'Create annotation key lists '
                           f'"{bed_key_conf}" and "{bw_key_conf}"')
        data_dir = getattr(self, 'data_dir_symlink')
        annot_key_conf = getattr(self, 'annot_key_conf')

        if annot_key_conf is None:
            create.create_annotation_key(bed_key_conf, data_dir, 'bed')
            create.create_annotation_key(bw_key_conf, data_dir, 'bw')
        else:
            create.split_annotation_key(bed_key_conf, bw_key_conf,
                                        annot_key_conf)

import argparse
from pathlib import Path

import cwas.core.configuration.create as create
import cwas.utils.log as log
from cwas.runnable import Runnable


class Configuration(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        work_dir = getattr(self, 'work_dir')
        self.data_dir_symlink = work_dir / 'annotation-data'
        self.gene_matrix_symlink = work_dir / 'gene_matrix.txt'
        self.bed_key_list = work_dir / 'annotation_key_bed.yaml'
        self.bw_key_list = work_dir / 'annotation_key_bw.yaml'
        self.bw_cutoff_list = work_dir / 'annotation_cutoff_bw.yaml'
        self.category_domain_list = work_dir / 'category_domain.yaml'
        self.redundant_category_table = work_dir / 'redundant_category.txt'

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
        self._create_data_dir_symlink()
        self._create_gene_matrix_symlink()
        self._create_annotation_key_list()
        self._create_bw_cutoff_list()
        self._create_category_info()
        self._set_env()

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
        bed_key_list = getattr(self, 'bed_key_list')
        bw_key_list = getattr(self, 'bw_key_list')
        log.print_progress(f'Create annotation key lists '
                           f'"{bed_key_list}" and "{bw_key_list}"')
        data_dir = getattr(self, 'data_dir_symlink')
        annot_key_conf = getattr(self, 'annot_key_conf')

        if annot_key_conf is None:
            create.create_annotation_key(bed_key_list, data_dir, 'bed')
            create.create_annotation_key(bw_key_list, data_dir, 'bw')
        else:
            create.split_annotation_key(bed_key_list, bw_key_list,
                                        annot_key_conf)

    def _create_bw_cutoff_list(self):
        bw_cutoff_list = getattr(self, 'bw_cutoff_list')
        bw_key_list = getattr(self, 'bw_key_list')
        bw_cutoff_conf = getattr(self, 'bw_cutoff_conf')  # User-defined
        log.print_progress(f'Create BigWig cufoff list "{bw_cutoff_list}"')
        create.create_bw_cutoff_list(bw_cutoff_list, bw_key_list,
                                     bw_cutoff_conf)

    def _create_category_info(self):
        """ Create a list of category domains and a redundant category table"""
        domain_list = getattr(self, 'category_domain_list')
        bed_key_list = getattr(self, 'bed_key_list')
        bw_key_list = getattr(self, 'bw_key_list')
        gene_matrix = getattr(self, 'gene_matrix')
        redundant_category_table = getattr(self, 'redundant_category_table')

        log.print_progress(f'Create a CWAS category domain list '
                           f'"{domain_list}"')
        create.create_category_domain_list(domain_list, bed_key_list,
                                           bw_key_list, gene_matrix)
        log.print_progress(f'Create a redundant category table '
                           f'"{redundant_category_table}"')
        create.create_redundant_category_table(redundant_category_table)

    def _set_env(self):
        log.print_progress('Set CWAS environment variables')
        cwas_env = getattr(self, 'env')
        cwas_env.set_env('CWAS_WORKSPACE', getattr(self, 'work_dir'))
        cwas_env.set_env('ANNOTATION_DATA', self.data_dir_symlink)
        cwas_env.set_env('GENE_MATRIX', self.gene_matrix_symlink)
        cwas_env.set_env('ANNOTATION_BED_KEY', self.bed_key_list)
        cwas_env.set_env('ANNOTATION_BW_KEY', self.bw_key_list)
        cwas_env.set_env('ANNOTATION_BW_CUTOFF', self.bw_cutoff_list)
        cwas_env.set_env('CATEGORY_DOMAIN', self.category_domain_list)
        cwas_env.set_env('REDUNDANT_CATEGORY', self.redundant_category_table)

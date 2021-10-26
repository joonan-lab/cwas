import argparse
import os

import cwas.utils.check as check
import cwas.utils.log as log
from cwas.runnable import Runnable


class Categorization(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._assign_config_to_attr('gene_mat_path', 'gene_matrix.txt')
        self._assign_config_to_attr('cat_conf_path', 'categories.yaml')
        self._assign_config_to_attr('rdd_cat_path', 'redundant_categories.yaml')

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True,
                            type=str,
                            help='Input VCF file from VEP')
        parser.add_argument('-o', '--outfile', dest='outfile_path',
                            required=False, type=str,
                            help='Path of the output',
                            default='cwas_categorization_result.txt')
        parser.add_argument('-p', '--num_proc', dest='num_proc', required=False,
                            type=int,
                            help='Number of worker processes for the '
                                 'categorization',
                            default=1)
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg('The input VCF file', args.in_vcf_path)
        log.print_arg('The output path', args.outfile_path)
        log.print_arg('No. worker processes for the categorization',
                      f'{args.num_proc: ,d}')

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        error.check_is_file(args.in_vcf_path)
        outfile_dir = os.path.dirname(args.outfile_path)
        error.check_is_dir(outfile_dir)
        error.check_num_proc(args.num_proc)

    def run(self):
        log.print_log('Notice', 'Not implemented yet.')

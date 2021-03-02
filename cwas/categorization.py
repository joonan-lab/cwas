import argparse
import os
from importlib.resources import path

import cwas.utils.error as error
import cwas.utils.log as log
from cwas.runnable import Runnable


class Categorization(Runnable):
    def _set_attr(self):
        with path('cwas.config', 'gene_matrix') as gene_mat_path:
            self.gene_mat_path = gene_mat_path
        with path('cwas.config', 'categories.yaml') as cat_conf_path:
            self.cat_conf_path = cat_conf_path
        with path('cwas.config', 'redundant_categories.yaml') as rdd_cat_path:
            self.rdd_cat_path = rdd_cat_path

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
        pass

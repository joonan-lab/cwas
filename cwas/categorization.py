import argparse
import multiprocessing as mp
import os

import cwas.utils.log as log
from cwas.runnable import Runnable


class Categorization(Runnable):
    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-i', '--infile', dest='in_vcf_path', required=True,
                            type=str,
                            help='Input VCF file from VEP')
        parser.add_argument('-o', '--outfile', dest='outfile_path',
                            required=False, type=str,
                            help='Path of the output',
                            default='cwas_cat_result.txt')
        parser.add_argument('-p', '--num_proc', dest='num_proc', required=False,
                            type=int,
                            help='Number of processes for this script',
                            default=1)

        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg('The input VCF file', args.in_vcf_path)
        log.print_arg('The output path', args.outfile_path)
        log.print_arg('No. processes for the categorization',
                      f'{args.num_proc: ,d}')

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        assert os.path.isfile(args.in_vcf_path), \
            f'The input VCF file "{args.in_vcf_path}" cannot be found.'
        outfile_dir = os.path.dirname(args.outfile_path)
        assert outfile_dir == '' or os.path.isdir(outfile_dir), \
            f'The outfile directory "{outfile_dir}" cannot be found.'
        assert 1 <= args.num_proc <= mp.cpu_count(), \
            f'Invalid number of processes "{args.num_proc:,d}". ' \
            f'It must be in the range [1, {mp.cpu_count()}].'

    def run(self):
        pass

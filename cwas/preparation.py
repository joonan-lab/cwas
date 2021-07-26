import argparse
import multiprocessing as mp
from cwas.runnable import Runnable
import cwas.utils.log as log


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
        min_num_proc = 1
        max_num_proc = mp.cpu_count()
        if args.num_proc < min_num_proc or args.num_proc > max_num_proc:
            raise ValueError(f'Wrong No. processes "{args.num_proc}" '
                             f'(range: {min_num_proc} ~ {max_num_proc})')

    def run(self):
        log.print_log('Notice', 'Not implemented yet.')

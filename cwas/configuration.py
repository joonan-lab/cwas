import argparse
from pathlib import Path

import cwas.utils.log as log
from cwas.runnable import Runnable


class Configuration(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description='Arguments for CWAS Configuration',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        default_work_dir = Path.home() / '.cwas'
        parser.add_argument('-d', '--work_dir', dest='work_dir', required=False,
                            type=Path, default=default_work_dir,
                            help='Path to your CWAS workspace'
                            )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
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
        try:
            self.args.work_dir.mkdir(parents=True, exist_ok=True)
        except NotADirectoryError:
            log.print_err('The input CWAS workspace path is invalid.')
            raise

        log.print_log('Notice', 'Not implemented yet.')

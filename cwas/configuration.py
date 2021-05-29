import argparse
import os
from pathlib import Path

import dotenv
import yaml

import cwas.config
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
        parser.add_argument('-d', '--annotation_data_dir', dest='data_dir',
                            required=True, type=Path,
                            help="Path to your annotation data directory"
                            )
        parser.add_argument('-c', '--bigwig_cutoff', dest='bw_cutoff_conf',
                            required=False, type=Path, default=None,
                            help="Path to a configuration file (.yaml) that "
                                 "specifies the annotation cutoff of "
                                 "each BigWig file")
        parser.add_argument('-k', '--annotation_key', dest='annot_key_conf',
                            required=False, type=Path, default=None,
                            help="Path to a configuration file (.yaml) that "
                                 "specifies the annotation key of each "
                                 "annotation data file")
        parser.add_argument('-w', '--workspace', dest='work_dir',
                            required=False,
                            type=Path, default=default_work_dir,
                            help='Path to your CWAS workspace'
                            )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg('CWAS workspace', args.work_dir)
        log.print_arg('Your annotation data directory', args.data_dir)

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
        work_dir = getattr(self, 'work_dir')
        data_dir = getattr(self, 'data_dir')
        data_dir_symlink = work_dir / 'annotation-data'
        try:
            log.print_progress(f'Create CWAS workspace "{work_dir}"')
            work_dir.mkdir(parents=True, exist_ok=True)
            log.print_progress(f'Create a symlink "{data_dir}" of your data '
                               f'directory')
            os.symlink(data_dir, data_dir_symlink, target_is_directory=True)
        except NotADirectoryError:
            log.print_err('The path to CWAS workspace is invalid.')
            raise
        except FileExistsError:
            log.print_warn(f'"{data_dir_symlink}" already exists so skip '
                           f'making symbolic link for your data directory.')
        data_dir = data_dir_symlink

        annot_key_conf = getattr(self, 'annot_key_conf')
        bed_key_conf = work_dir / 'annotation_key_bed.yaml'
        bw_key_conf = work_dir / 'annotation_key_bw.yaml'
        if annot_key_conf is None:
            log.print_progress('Create a annotation key list (.yaml)')
            bed_filenames = [str(bed_filepath)
                             for bed_filepath in data_dir.glob('*.bed.gz')]
            bed_key_dict = {bed_filename: bed_filename[:-7].replace('.', '_')
                            for bed_filename in bed_filenames}
            bw_filenames = [str(bw_filepath)
                            for bw_filepath in data_dir.glob('*.bw')]
            bw_key_dict = {bw_filename: bw_filename[:-3].replace('.', '_')
                           for bw_filename in bw_filenames}
            with bed_key_conf.open('w') as bed_key_f:
                yaml.dump(bed_key_dict, bed_key_f)
            with bw_key_conf.open('w') as bw_key_f:
                yaml.dump(bw_key_dict, bw_key_f)

        log.print_progress('Create a dotenv')
        cwas_config_dir = Path(cwas.config.__file__).parent
        env_path = cwas_config_dir.resolve() / '.env'
        env_path.touch()
        dotenv.set_key(env_path, 'CWAS_WORKSPACE', str(work_dir))

        log.print_log('Notice', 'Not implemented yet.')

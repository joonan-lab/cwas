"""
Functions to make and execute system commands
"""
import os
import sys
import cwas.utils.log as log


def execute(cmd: str, raise_err: bool = False):
    log.print_log('CMD', cmd, True)
    exit_val = os.system(cmd)

    if exit_val != 0:
        # This command has failed.
        msg = f'This CMD has failed with a exit value {exit_val}.'
        if raise_err:
            log.print_err(msg)
            sys.exit()
        else:
            log.print_warn(msg)


def bgzip_tabix(in_file_path: str, force_overwrite: int = 0):
    """ Block compression (bgzip) and make an index (tabix)

    :param in_file_path: Path of a TAB-delimited genome position files
    :param force_overwrite: This parameter will be considered as boolean.
                            If 1, the bgzip and tabix command will be executed
                            regardless of existence of the result files.
    """
    gz_path = in_file_path + '.gz'

    if not force_overwrite and os.path.isfile(gz_path):
        log.print_progress('"{in_file_path}" has already been bgzipped so '
                           'skip the bgzip step.')
    else:
        log.print_progress(f'Execute bgzip and tabix for "{in_file_path}".')
        cmd = f'bgzip {in_file_path};'
        cmd += f'tabix {gz_path};'
        execute(cmd)

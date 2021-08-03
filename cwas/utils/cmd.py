"""
Functions to make and execute system commands
"""
import shutil
import subprocess

import cwas.utils.log as log


def execute(args: list, raise_err: bool = True) -> subprocess.CompletedProcess:
    cmd = ' '.join(args)
    log.print_log('CMD', cmd, True)

    try:
        output = subprocess.run(args, check=raise_err)
        if output.returncode != 0:
            msg = f'This CMD has failed with exit status {output.returncode}.'
            log.print_warn(msg)
        return output
    except subprocess.CalledProcessError:
        log.print_err('Your command has failed.')
        raise


def execute_bin(
    bin_name: str,
    args: list,
    raise_err: bool = True
) -> subprocess.CompletedProcess:
    bin_path = shutil.which(bin_name)

    if bin_path is None:
        raise FileNotFoundError(
            f'"{bin_name}" cannot be found.')

    cmd = [bin_path] + args
    return execute(cmd, raise_err)

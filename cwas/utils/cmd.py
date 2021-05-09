"""
Functions to make and execute system commands
"""
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

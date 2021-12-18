"""
Utils to execute system command line.
"""
import shutil
import subprocess
from pathlib import Path

import cwas.utils.log as log
from cwas.utils.check import check_is_file


class CmdExecutor:
    def __init__(self, bin_name: str, args: list = []):
        self._bin_name = bin_name
        self._args = list(args)

    @property
    def cmd(self) -> list:
        return [self._bin_name, *self._args]

    def execute(self) -> int:
        log.print_log("CMD", " ".join(self.cmd), True)
        output = subprocess.run(self.cmd, check=False)

        if output.returncode:
            log.print_warn(
                f"Command {self.cmd} has failed"
                f"with exit status {output.returncode}."
            )
        return output.returncode

    def execute_raising_err(self) -> int:
        log.print_log("CMD", " ".join(self.cmd), True)
        output = subprocess.run(self.cmd, check=True)
        return output.returncode


def execute(args: list, raise_err: bool = True) -> subprocess.CompletedProcess:
    cmd = " ".join(args)
    log.print_log("CMD", cmd, True)

    try:
        output = subprocess.run(args, check=raise_err)
        if output.returncode != 0:
            msg = f"This CMD has failed with exit status {output.returncode}."
            log.print_warn(msg)
        return output
    except subprocess.CalledProcessError:
        log.print_err("Your command has failed.")
        raise


def execute_bin(
    bin_name: str, args: list, raise_err: bool = True
) -> subprocess.CompletedProcess:
    bin_path = shutil.which(bin_name)

    if bin_path is None:
        raise FileNotFoundError(f'The binary "{bin_name}" cannot be found.')

    cmd = [bin_path] + args
    return execute(cmd, raise_err)


def compress_bed_file(bed_file_path: Path) -> Path:
    """Compress the BED file using bgzip"""
    gz_path = Path(str(bed_file_path) + ".gz")
    if gz_path.exists():
        log.print_warn(
            f'The compressed BED file "{bed_file_path}" already exists. '
            f"Skip compressing."
        )
    else:
        check_is_file(bed_file_path)

        try:
            execute_bin("bgzip", [str(bed_file_path)])
        except subprocess.CalledProcessError:
            log.print_err(
                f'Failed to compress your BED file "{bed_file_path}".'
            )
            raise
        except FileNotFoundError:
            log.print_err('"bgzip" is not installed in your environment.')
            raise

    return gz_path


def index_bed_file(comp_bed_path: Path) -> Path:
    """Index the compressed BED file"""
    idx_path = Path(str(comp_bed_path) + ".tbi")
    if idx_path.exists():
        log.print_warn(
            f'The BED file index "{comp_bed_path}" already exists. '
            f"Skip indexing."
        )
    else:
        check_is_file(comp_bed_path)

        try:
            execute_bin("tabix", [str(comp_bed_path)])
        except subprocess.CalledProcessError:
            log.print_err(
                f'Failed to index your compressed BED file "{comp_bed_path}".'
            )
            raise
        except FileNotFoundError:
            log.print_err('"tabix" is not installed in your environment.')
            raise

    return idx_path

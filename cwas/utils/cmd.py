"""
Utils to execute system command line.
"""
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import cwas.utils.log as log
from cwas.utils.check import check_is_file


class CmdExecutor:
    def __init__(self, bin: str, args: list = []):
        def resolve_bin(bin: str) -> Optional[str]:
            if Path(bin).resolve().is_file():
                return str(Path(bin).resolve())
            return shutil.which(bin)

        self._bin_path = resolve_bin(bin)
        if self._bin_path is None:
            raise FileNotFoundError(f"Failed to find the binary '{bin}'")
        self._args = list(args)

    @property
    def bin_path(self):
        return self._bin_path

    @property
    def cmd(self) -> list:
        return [self._bin_path, *self._args]

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
            CmdExecutor("bgzip", [str(bed_file_path)]).execute_raising_err()
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
            CmdExecutor("tabix", [str(comp_bed_path)]).execute_raising_err()
        except subprocess.CalledProcessError:
            log.print_err(
                f'Failed to index your compressed BED file "{comp_bed_path}".'
            )
            raise
        except FileNotFoundError:
            log.print_err('"tabix" is not installed in your environment.')
            raise

    return idx_path

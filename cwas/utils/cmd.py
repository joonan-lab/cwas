"""
Utils to execute system command line.
"""
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import cwas.utils.log as log


class CmdExecutor:
    def __init__(self, bin: str, args: list = [], multi_input: Optional[str] = None, shell: bool = False):
        def resolve_bin(bin: str) -> Optional[str]:
            if Path(bin).resolve().is_file():
                return str(Path(bin).resolve())
            return shutil.which(bin)

        self._bin_path = resolve_bin(bin)
        if self._bin_path is None:
            raise FileNotFoundError(f"Failed to find the binary '{bin}'")
        self._args = list(args)
        self._multi_input = multi_input if multi_input is not None else None
        self._shell = shell

    @property
    def bin_path(self):
        return self._bin_path

    @property
    def cmd(self) -> list:
        if self._multi_input is None:
            cmd_args_list = [self._bin_path, *self._args]
        else:
            cmd_args_list = [self._multi_input, self._bin_path, *self._args]
        return cmd_args_list

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
        if self._shell == True:
            log.print_log("CMD", " ".join(self.cmd), True)
            output = subprocess.run(" ".join(self.cmd), check=True, shell = True)            
        else:
            log.print_log("CMD", " ".join(self.cmd), True)
            output = subprocess.run(self.cmd, check=True, shell = False)
        return output.returncode

def compress_using_bgzip(
    input_file_path: Path, force_overwrite: bool = False
) -> Path:
    result = Path(str(input_file_path) + ".gz")
    if not force_overwrite and result.exists():
        log.print_warn(
            f'The compressed file "{result}" already exists. '
            f"Skip compressing."
        )
    else:
        args = [str(input_file_path)]
        if force_overwrite:
            args.append("-f")
        CmdExecutor("bgzip", args).execute_raising_err()

    return result


def index_using_tabix(
    compressed_file_path: Path, force_overwrite: bool = False
) -> Path:
    result = Path(str(compressed_file_path) + ".tbi")
    if not force_overwrite and result.exists():
        log.print_warn(
            f'The file index "{result}" already exists. Skip indexing.'
        )
    else:
        CmdExecutor("tabix", [str(compressed_file_path)]).execute_raising_err()

    return result

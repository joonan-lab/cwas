"""
Utilities for the CWAS preparation step
"""
import subprocess
from pathlib import Path

import cwas.utils.log as log
from cwas.utils.cmd import execute_bin
from cwas.utils.error import check_is_file


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

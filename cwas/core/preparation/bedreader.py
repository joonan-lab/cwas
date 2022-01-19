"""
A container class that contains iterator protocol to read a BED file
"""
from pathlib import Path

import pysam
from cwas.utils.log import print_warn


class BedReader:
    def __init__(self, bed_path: Path):
        self.bed_path = bed_path
        self.contig = None
        self.start = None
        self.stop = None

        if not self.bed_path.exists():
            raise OSError(
                f"Failed to initialize a BedReader instance. "
                f"'{self.bed_path}' does not exists."
            )

    def __iter__(self):
        try:
            with pysam.TabixFile(str(self.bed_path)) as bed_file:
                for fields in bed_file.fetch(
                    self.contig, self.start, self.stop, parser=pysam.asTuple()
                ):
                    contig, start, stop, *others = tuple(fields)
                    yield (contig, int(start), int(stop), *others)
        except ValueError:
            print_warn(
                f"The contig '{self.contig}' does not exist "
                f"in the bed file '{self.bed_path}'."
            )

    def set_contig(self, contig: str):
        self.contig = contig

import argparse
from pathlib import Path


def categorization() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of CWAS categorization step",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Annotated VCF file",
    )
    result.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the categorization",
        default=1,
    )
    return result

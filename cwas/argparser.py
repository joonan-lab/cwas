import argparse
from pathlib import Path


def start() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments for Initializing a CWAS workspace",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = Path.home() / ".cwas"
    result.add_argument(
        "-w",
        "--workspace",
        dest="workspace",
        required=False,
        type=Path,
        default=default_workspace,
        help="Path to your CWAS workspace directory",
    )
    return result


def preparation() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments for Annotation Data Preparation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        default=1,
        help="Max No. processes for this step",
    )
    result.add_argument(
        "-f",
        "--force_overwrite",
        dest="force_overwrite",
        action="store_const",
        const=1,
        default=0,
        help="Force to overwrite the result",
    )
    return result


def annotation() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of CWAS annotation step",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-v",
        "--vcf_file",
        dest="vcf_path",
        required=True,
        type=Path,
        help="Target VCF file",
    )
    return result


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


def binomial_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of Burden Tests",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file",
    )
    result.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=True,
        type=Path,
        help="File listing information of your samples",
    )
    result.add_argument(
        "-a",
        "--adjustment_factor",
        dest="adj_factor_path",
        required=False,
        default=None,
        type=Path,
        help="File listing adjustment factors of each sample",
    )
    return result

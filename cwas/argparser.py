import argparse
from pathlib import Path

import dotenv


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

def configuration() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments for CWAS Configuration",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-d",
        "--annotation_data_dir",
        dest="data_dir",
        required=False,
        type=Path,
        help="Path to your annotation data directory",
    )
    result.add_argument(
        "-m",
        "--gene_matrix",
        dest="gene_matrix",
        required=False,
        type=Path,
        help="Path to your gene matrix",
    )
    result.add_argument(
        "-a",
        "--annotation_key_config",
        dest="annot_key_conf",
        required=False,
        type=Path,
        help="Path to a configuration file (.yaml) that "
        "specifies the annotation key of each "
        "annotation data file",
    )
    result.add_argument(
        "-v",
        "--vep",
        dest="vep",
        required=False,
        type=Path,
        help="Path to Variant Effect Predictor (VEP)",
    )
    result.add_argument(
        "-vrd",
        "--vep_resource_dir",
        dest="vep_resource_dir",
        required=False,
        type=Path,
        help="Path to your VEP resource directory",
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
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-v",
        "--vcf_file",
        dest="vcf_path",
        required=True,
        type=Path,
        help="Target VCF file",
    )
    result.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the categorization (default: 1)",
        default=1,
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved",
    )
    return result


def categorization() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of CWAS categorization step",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Annotated VCF file",
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved",
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
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file",
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved",
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
    result.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        required=False,
        action="store_true",
        help="Use the number of samples with variants in each category for burden test instead of the number of variants",
    )
    return result


def permutation_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of Burden Tests",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file (gzipped)",
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved",
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
    result.add_argument(
        "-n",
        "--num_perm",
        dest="num_perm",
        default=10000,
        type=int,
        help="The number of label-swapping permutations",
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
    result.add_argument(
        "-b",
        "--burden_shift",
        dest="burden_shift",
        required=False,
        action="store_true",
        help="Generate a file of binomial p-values for each burden-shifted data",
        )
    result.add_argument(
        "-rr",
        "--perm_rr",
        dest="save_perm_rr",
        required=False,
        action="store_true",
        help="Generate a file of relative risks (RRs) for each burden-shifted data",
    )
    result.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        required=False,
        action="store_true",
        help="Use the number of samples with variants in each category for burden test instead of the number of variants",
    )
    return result

def extract_variant() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of Variant Extraction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Annotated VCF file",
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved",
    )
    result.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output file (i.e., output.<tag>.extracted_variants.txt.gz)",
    )
    result.add_argument(
        "-c",
        "--category_set_path",
        dest="category_set_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file containing categories for extracting variants",
    )
    result.add_argument(
        "-ai",
        "--annotation_info",
        dest="annotation_info",
        required=False,
        default=False,
        action="store_true",
        help="Save with annotation information attached (such as gene list, functional annotations, etc)",
    )    
    return result

def simulation() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of random variant generation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        '-i', '--in_vcf',
        dest='in_vcf_path',
        required=True,
        type=Path,
        help='Input VCF file which is referred to generate random mutations'
    )
    result.add_argument(
        '-s',
        '--sample_info',
        dest='sample_info_path',
        required=True,
        type=Path,
        help='File listing sample IDs with their families and sample_types (case or ctrl)'
    )
    result.add_argument(
        '-o',
        '--out_dir',
        dest='out_dir',
        required=False,
        type=Path,
        help='Directory of outputs that lists random mutations. '
        'The number of outputs will be the same with the number of simulations. '
        '(default: $CWAS_WORKSPACE/random-mutation)'
    )
    result.add_argument(
        '-t',
        '--out_tag',
        dest='out_tag',
        required=False,
        type=str,
        help='Prefix of output files. Each output file name will start with this tag.',
        default='rand_mut'
    )
    result.add_argument(
        '-n',
        '--num_sim',
        dest='num_sim',
        required=False,
        type=int,
        help='Number of simulations to generate random mutations',
        default=1
    )
    result.add_argument(
        '-p',
        '--num_proc',
        dest='num_proc',
        required=False,
        type=int,
        help='Number of processes for this script (only necessary for split VCF files)',
        default=1
    )
    result.add_argument(
        "-r",
        "--resume",
        dest="resume",
        required=False,
        default=False,
        action="store_true",
        help="Resume the simulation from the last step. Assume some generated output files are not truncated."
    )
    return result

def multiprocessing() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of CWAS multiprocessing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    def add_common_args(subparser: argparse.ArgumentParser):
        """ Add common arguments to the subparser """
        subparser.add_argument('-cs', '--cwas_step', dest='cwas_step', required=True, type=str, default='all',
                            help='CWAS process for multiprocessing', choices=['all', 'annotation', 'categorizatoin', 'binomial_test'])
        subparser.add_argument('-i', '--input_dir', dest='in_dir_path', required=True, type=Path,
                            help='Path to the input directory')
        subparser.add_argument(
            "-o_dir",
            "--output_directory",
            dest="output_dir_path",
            required=False,
            type=Path,
            help="Directory where output file will be saved (default : $CWAS_WORKSPACE/random-mutation_{{args.step}})",
        )
        subparser.add_argument('-mp', '--num_mp', dest='num_mp', required=False, type=int,
                            help='Number of concurrent processes for CWAS multiprocessing', default=1)
        subparser.add_argument(
            "-r", "--resume",
            dest="resume", required=False, default=False, action="store_true",
            help="Resume the simulation from the last step. Assume some generated output files are not truncated.",
        )
        
    subparsers = result.add_subparsers(description='A name of a step of CWAS {annotation, categorization, binomial_test}',
                                       dest='step')
    parser_annot = subparsers.add_parser(
        'annotation',
        description='Multiprocessing variant annotation in CWAS',
        help='Multiprocessing variant annotation in CWAS (arg "annotation -h" for usage)'
    )
    parser_annot.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the annotation (default: 1)",
        default=1,
    )
    add_common_args(parser_annot)
        
    parser_cat = subparsers.add_parser(
        'categorization',
        description='Multiprocessing variant categorization in CWAS',
        help='Multiprocessing variant categorization in CWAS (arg "categorization -h" for usage)'
    )
    parser_cat.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the categorization (default: 1)",
        default=1,
    )
    add_common_args(parser_cat)

    parser_burden = subparsers.add_parser(
        'binomial_test',
        description='Multiprocessing burden binomial tests in CWAS',
        help='Multiprocessing burden binomial tests in CWAS (arg "binomial_test -h" for usage)'
    )
    parser_burden.add_argument('-s', '--sample_file', dest='sample_file_path', required=True, type=str,
                               help='File listing sample IDs with their families and sample_types (case or ctrl)')
    parser_burden.add_argument('-a', '--adj_file', dest='adj_file_path', required=False, type=str,
                               help='File that contains adjustment factors for No. DNVs of each sample', default='')
    parser_burden.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        required=False,
        default=False,
        action="store_true",
        help="Use the number of samples with variants in each category for burden test instead of the number of variants",
    )
    add_common_args(parser_burden)
    
    return result
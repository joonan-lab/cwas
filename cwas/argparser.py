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
    parser.add_argument(
        "-vmdb",
        "--vep_mis_db",
        dest="vep_mis_db",
        required=False,
        type=Path,
        help="Path to your VCF file (database for predicting damaging missense mutations)",
    )
    parser.add_argument(
        "-vmk",
        "--vep_mis_info_key",
        dest="vep_mis_info_key",
        required=False,
        type=Path,
        help="VCF info field key name in  database",
    )
    parser.add_argument(
        "-vmt",
        "--vep_mis_thres",
        dest="vep_mis_thres",
        required=False,
        type=float,
        help="Threshold of score from database for predicting damaging missense mutations",
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
        '-p',
        '--num_proc',
        dest='num_proc',
        required=False,
        default=1,
        type=int,
        help="Number of processes for the annotation (default: 1)",
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
    result.add_argument(
        "-m",
        "--matrix",
        dest="generate_matrix",
        required=False,
        action="store_true",
        help="Generate a correlation matrix and a matrix with intersected number of variants bewteen categories",
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
        help="Number of worker processes for the permutation",
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
        help='Number of simulations to generate random variants',
        default=1
    )
    result.add_argument(
        '-p',
        '--num_proc',
        dest='num_proc',
        required=False,
        type=int,
        help='Number of processes for random variant generation (only necessary for split VCF files)',
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


def concat_zscore() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of Z-score concatenation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i_dir",
        "--input_directory",
        dest="input_dir_path",
        required=True,
        type=Path,
        help="Directory where burden test results are saved",
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
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output file (i.e., output.<tag>.zscores.txt.gz)",
    )
    result.add_argument(
        "-c",
        "--category_set_path",
        dest="category_set_path",
        required=True,
        default=None,
        type=Path,
        help="Path to a text file containing categories for extracting variants",
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
    return result


def effective_num_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of Effective Number of Test Calculation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Path to the input file",
    )
    result.add_argument(
        "-if",
        "--input_format",
        dest="input_format",
        required=False,
        default = 'corr',
        choices = ['corr', 'inter', 'zscores'],
        type=str,
        help="Input format. If not specified, corr will be used.\nAvailable options:\n\tcorr: a correlation matrix\n\tinter: a matrix with intersected number of variants between categories\n\tzscores: concatenated z-scores",
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
        '-n',
        '--num_sim',
        dest='num_sim',
        required=False,
        type=int,
        help='Number of eigen values to use',
        default=10000
    )
    result.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=False,
        type=Path,
        help="File listing information of your samples. Required only when input format is set to 'inter'",
    )
    result.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output files (e.g., corr_mat_<tag>.pickle)",
    )
    result.add_argument(
        "-c",
        "--category_set_path",
        dest="category_set_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file containing categories for eigen decomposition. If not specified, all of the categories in the z-score file will be used.",
    )
    result.add_argument(
        "-ef",
        "--eff_num_test",
        dest="eff_num_test",
        required=False,
        action="store_true",
        help="Calculate the effective number of tests",
    )
    return result

def dawn() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Arguments of CWAS DAWN analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    result.add_argument(
        "-i",
        "--input_directory",
        dest="input_dir_path",
        required=True,
        type=Path,
        help="Directory where input files stored. This directory must include three required input files.\n 1. Eigen vectors file (*eig_vecs*.txt.gz)\n 2. Category correlation matrix file (*cor_mat*.pickle)\n 3. Permutation test file (*permutation_test*.txt.gz)",
    )
    result.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        default=1,
        type=int,
        help="Number of worker processes for the DAWN.",
    )
    result.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved.",
    )
    result.add_argument(
        "-r",
        "--range",
        dest="k_range",
        required=False,
        default='2,100',
        type=str,
        help="Range from start and end to find K for k-means clustering, and start must be above 1. Start and end should be comma-separated. -r, and -k are mutually exclusive.\nIf you want to use range (-r) to find optimal K, -k must be None.",
    )
    result.add_argument(
        "-k",
        "--k_val",
        dest="k_val",
        required=False,
        type=int,
        default=None,
        help="K for k-means clustering. -r, and -k are mutually exclusive. If -k is given, the value of -r is ignored.",
    )
    result.add_argument(
        "-s",
        "--seed",
        dest="seed",
        required=False,
        default=1,
        type=int,
        help="Seed value for t-SNE.",
    )
    result.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=True,
        type=str,
        help="Tag used for the name of output files (e.g. intergenic, coding etc.).",
    )
    result.add_argument(
        "-c_set",
        "--categroy_file",
        dest="category_set_file",
        required=True,
        type=Path,
        help="File path of category set file used for DAWN analysis.",
    )
    result.add_argument(
        "-c_count",
        "--cat_count",
        dest="category_count_file",
        required=True,
        type=Path,
        help="File path of category counts file resulted from burden test (for each variant) or sign test (for each sample).",
    )
    result.add_argument(
        "-D",
        "--count_threshold",
        dest="count_threshold",
        required=False,
        type=int,
        default=20,
        help="The treshold of DNM counts. The least amount of variants a category should have.",
    )
    result.add_argument(
        "-C",
        "--corr_threshold",
        dest="corr_threshold",
        required=False,
        type=float,
        default=0.12,
        help="The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster.",
    )
    result.add_argument(
        "-S",
        "--size_threshold",
        dest="size_threshold",
        required=False,
        type=int,
        default=2,
        help="The threshold of the number of categories per cluster. The least amount of categories a cluster should have.",
    )
            
    return result
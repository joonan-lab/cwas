import argparse
from pathlib import Path

import dotenv


def start() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Initializing a CWAS workspace",
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
        description="CWAS Configuration",
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
    result.add_argument(
        "-vmdb",
        "--vep_mis_db",
        dest="vep_mis_db",
        required=False,
        type=Path,
        help="Path to your VCF file (database for predicting damaging missense mutations)",
    )
    result.add_argument(
        "-vmk",
        "--vep_mis_info_key",
        dest="vep_mis_info_key",
        required=False,
        type=Path,
        help="VCF info field key name in  database",
    )
    result.add_argument(
        "-vmt",
        "--vep_mis_thres",
        dest="vep_mis_thres",
        required=False,
        type=float,
        help="Threshold of score from database for predicting damaging missense mutations",
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

def preparation() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas preparation",
        description="Preparation of Annotation Data for CWAS",
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
        prog="cwas annotation",
        description="CWAS annotation step",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-v",
        "--vcf_file",
        dest="vcf_path",
        required=True,
        type=Path,
        help="Target VCF file",
    )
    optional.add_argument(
        '-p',
        '--num_proc',
        dest='num_proc',
        required=False,
        default=1,
        type=int,
        help="Number of processes for the annotation (default: 1)",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {})".format(default_workspace),
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    return result


def categorization() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas categorization",
        description="CWAS categorization step",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Annotated VCF file",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {})".format(default_workspace),
    )
    optional.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the categorization (default: 1)",
        default=1,
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    return result


def binomial_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas binomial_test",
        description="CWAS burden tests step - binomial test",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file (*.zarr) resulted from categorization step.",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved. (default: {})".format(default_workspace),
    )
    required.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=True,
        type=Path,
        help="File listing information of your samples.",
    )
    optional.add_argument(
        "-a",
        "--adjustment_factor",
        dest="adj_factor_path",
        required=False,
        default=None,
        type=Path,
        help="File listing adjustment factors of each sample. The file is required to use the adjusted values in the binomial test. (default: None)",
    )
    optional.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        required=False,
        action="store_true",
        help="Use the number of samples with variants in each category for burden test instead of the number of variants.",
    )
    optional.add_argument(
        '-t',
        '--tag',
        dest='tag',
        default=None,
        type=str,
        required=False,
        help="Tags of category queried for highlighting points on the volcano plot. If you use multiple tags, concatenate by ','. (e.g. CRE,CHD8) (default: None)"
    )
    optional.add_argument(
        "-num_ef",
        "--num_effective_test",
        dest="eff_test",
        required=False,
        type=int,
        default=False,
        help="Number of effective tests (default: False)",
    )
    optional.add_argument(
        "-ms",
        "--marker_size",
        dest="marker_size",
        required=False,
        type=float,
        default=15,
        help="Maker size of the volcano plot resulted from the binomial test. (unit: pt) (default: 15)",
    )
    optional.add_argument(
        "-fs",
        "--font_size",
        dest='font_size',
        required=False,
        type=float,
        default=15,
        help="Font size of the volcano plot resulted from the binomial test. (unit: pt) (default: 15)",
    )
    optional.add_argument(
        '-ps',
        '--plot_size',
        dest="plot_size",
        required=False,
        type=float,
        default=7,
        help="Plot size of the volcano plot resulted from the binomial test, width and height are the same. (unit: inch) (default: 7)"
    )
    optional.add_argument(
        "-pt",
        "--plot_title",
        dest="plot_title",
        required=False,
        type=str,
        default="Binomial test result",
        help="Title of volcano plot of binomial test result (default: Binomial test result)."
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result


def permutation_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas permutation_test",
        description="CWAS burden tests step - permutation test",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file (gizpped) resulted from categorization step.",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved. (default: {})".format(default_workspace),
    )
    required.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=True,
        type=Path,
        help="File listing information of your samples",
    )
    optional.add_argument(
        "-a",
        "--adjustment_factor",
        dest="adj_factor_path",
        required=False,
        default=None,
        type=Path,
        help="File listing adjustment factors of each sample. The file is required to use the adjusted values in the binomial test. (default: None)",
    )
    optional.add_argument(
        "-n",
        "--num_perm",
        dest="num_perm",
        default=10000,
        type=int,
        help="The number of label-swapping permutations. (default: 10000)",
    )
    optional.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the permutation. (default: 1)",
        default=1,
    )
    optional.add_argument(
        "-b",
        "--burden_shift",
        dest="burden_shift",
        required=False,
        action="store_true",
        help="Generate a file of binomial p-values for each burden-shifted data",
        )
    #optional.add_argument(
    #    "-rr",
    #    "--perm_rr",
    #    dest="save_perm_rr",
    #    required=False,
    #    action="store_true",
    #    help="Generate a file of relative risks (RRs) for each burden-shifted data",
    #)
    optional.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        required=False,
        action="store_true",
        help="Use the number of samples with variants in each category for burden test instead of the number of variants",
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result

def extract_variant() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas extract_variant",
        description="Extract interesting variants from annotated VCF file",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Annotated VCF file",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {})".format(default_workspace),
    )
    optional.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output file (i.e., output.<tag>.extracted_variants.txt.gz) (default: None)",
    )
    result.add_argument(
        '-c_set',
        '--category_set',
        dest="category_set_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file containing categories for extracting variants (default: None)",
    )
    result.add_argument(
        "-ai",
        "--annotation_info",
        dest="annotation_info",
        required=False,
        default=False,
        action="store_true",
        help="Save with annotation information attached (such as gene list, functional annotations, etc) (default: False)",
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result

def effective_num_test() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas effective_num_test",
        description="Calculation of Effective Number of Tests and Eigen Decomposition",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Path to the input file, either correlation matrix or intersection matrix resulting from categorization.",
    )
    required.add_argument(
        "-c_count",
        "--cat_count",
        dest="category_count_file",
        required=True,
        type=Path,
        help="File path of category counts file resulted from burden test (for each variant) or sign test (for each sample).",
    )
    optional.add_argument(
        "-if",
        "--input_format",
        dest="input_format",
        required=False,
        default = 'corr',
        choices = ['corr', 'inter'],
        type=str,
        help="Input format. If not specified, 'corr' will be used.\n"\
             "Available options:\n"\
             "* corr: a correlation matrix\n"\
             "* inter: a matrix with intersected number of variants between categories",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved. (default: {})".format(default_workspace),
    )
    optional.add_argument(
        '-n',
        '--num_eig',
        dest='num_eig',
        required=False,
        type=int,
        help='Number of eigen values to use. (default: 10000)',
        default=10000
    )
    optional.add_argument(
        "-thr",
        "--threshold",
        dest="count_thres",
        required=False,
        default=None,
        type=int,
        help="The number of variants (or samples) to filter categories (counts ≥ threshold)",
    )
    optional.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=False,
        type=Path,
        help="File listing information of your samples to calculate threshold of the number of variants (samples) of given categories. Required when '-thr' is not given.",
    )
    optional.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output files (e.g., <tag>.correlation_matrix.pkl). (default: None)",
    )
    optional.add_argument(
        "-c_set",
        "--category_set",
        dest="category_set_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file containing categories for eigen decomposition. If not specified, all of the categories will be used. (default: None)",
    )
    optional.add_argument(
        '-c_info',
        '--category_info',
        dest="category_info_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file with category information (*.category_info.txt).",
    )
    optional.add_argument(
        '-d',
        '--domain_list',
        dest="domain_list",
        required=False,
        default='all',
        type=str,
        help="Domain list to filter categories based on GENCODE domain. If 'run_all' is given, all available options will be tested (default: all).\n"\
             "Available options: run_all,all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA",
    )
    optional.add_argument(
        "-ef",
        "--eff_num_test",
        dest="eff_num_test",
        required=False,
        action="store_true",
        help="Calculate and output the effective number of tests. Only eigenvalues are generated with this option.",
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result


def burden_shift() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas burden_shift",
        description="CWAS Burden Shift Test step",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="input_path",
        required=True,
        type=Path,
        help="Path to the input file which is the result of burden test from binomial test (*.burden_test.txt)",
    )
    required.add_argument(
        '-b',
        '--burden_res',
        dest='burden_res',
        required=True,
        type=Path,
        help='Path to the result of burden shift from permutation test (*.binom_pvals.txt.gz)',
    )
    required.add_argument(
        '-c_info',
        '--category_info',
        dest='cat_set_file',
        required=True,
        type=Path,
        help='Path to a text file with category information (*.category_info.txt).',
    )
    required.add_argument(
        '-c_count',
        '--cat_count',
        dest='cat_count_file',
        required=True,
        type=Path,
        help='Path of the categories counts file from permutation test (*.category_counts.txt).',
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {})".format(default_workspace),
    )
    optional.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of the output files (default: None).",
    )
    optional.add_argument(
        "-c_cutoff",
        "--count_cutoff",
        dest="count_cutoff",
        required=False,
        default=7,
        type=int,
        help="The number of cutoff for category counts. It must be positive value (default: 7).",
    )
    optional.add_argument(
        "-pval",
        "--pval",
        dest="pval",
        required=False,
        default=0.05,
        type=float,
        help="P-value of threshold (default: 0.05).",
    )
    optional.add_argument(
        "-c_set",
        "--category_set",
        dest="cat_set_list",
        required=False,
        type=Path,
        help="Path of the list of interest category sets for the main output plot. In the file, one line stores one category set name and, do not include header.\nIf the category set name is combination of two or more domains, it must be separated by &. \nIf no user input is entered, the plot outputs for top N category sets."
    )
    optional.add_argument(
        "-N",
        "--n_cat_sets",
        dest='n_cat_sets',
        required=False,
        default=10,
        type=int,
        help="The number of the category sets contained in the main output plot. Top N category sets will be contain in main output plot (default: 10)."   
    )
    optional.add_argument(
        "-pt",
        "--plot_title",
        dest="plot_title",
        required=False,
        type=str,
        default='Burdenshift: Overrepresented terms',
        help="Title of summarized plot of burden shift result (default: Burdenshift: Overrepresented terms)."
    )
    optional.add_argument(
        "-fs",
        "--fontsize",
        dest='fontsize',
        required=False,
        default=10,
        type=int,
        help="Font size of final main output plot (default: 10)."   
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result

def risk_score() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas risk_score",
        description="CWAS risk score analysis",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="categorization_result_path",
        required=True,
        type=Path,
        help="The path of the categorization result file (*.zarr)",
    )
    required.add_argument(
        "-s",
        "--sample_info",
        dest="sample_info_path",
        required=True,
        type=Path,
        help="File listing sample IDs with their families and sample_types (case or ctrl).\nIf test categorization result is not available, 'SET' column is used for dividing training and test set."
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {})".format(default_workspace),
    )
    optional.add_argument(
        "-a",
        "--adjustment_factor",
        dest="adj_factor_path",
        required=False,
        default=None,
        type=Path,
        help="File listing adjustment factors of each sample. The file is required to use the adjusted values in the binomial test. (default: None)",
    )
    optional.add_argument(
        '-c_info',
        '--category_info',
        dest="category_set_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file with category information (*.category_info.txt).",
    )
    optional.add_argument(
        '-d',
        '--domain_list',
        dest="domain_list",
        required=False,
        default='all',
        type=str,
        help="Domain list to filter categories based on GENCODE domain. If 'run_all' is given, all available options will be tested (default: all).\n"\
             "Available options: run_all,all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA",
    )
    optional.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default='all',
        type=str,
        help="Tag used for the name of output files",
    )
    optional.add_argument(
        "-loop",
        "--do_loop",
        dest="do_loop",
        action="store_true",
        help="Use each annotation from functional annotation to calculate risk score.",
    )
    optional.add_argument(
        "-n_one",
        "--n_of_one_leave",
        dest="n_of_one_leave",
        action="store_true",
        help="Calculate risk score while excluding one annotation from functional annotation. This option is not used when the '--do_loop' flag is enabled.",
    )
    optional.add_argument(
        "-u",
        "--use_n_carrier",
        dest="use_n_carrier",
        action="store_true",
        help="Use the number of samples with variants in each category for calculating R2 instead of the number of variants.",
    )
    optional.add_argument(
        "-thr",
        "--threshold",
        dest="ctrl_thres",
        required=False,
        default=3,
        type=int,
        help="The number of variants in controls (or the number of control carriers) used to select rare categories (defulat: 3).",
    )
    optional.add_argument(
        '-tf',
        '--train_set_fraction',
        dest='train_set_f',
        required=False,
        type=float,
        default=0.7,
        help='Fraction of the training set (default: 0.7)')
    optional.add_argument(
        '-n_reg',
        '--num_regression',
        dest='num_reg',
        required=False,
        type=int,
        default=10,
        help='No. regression trials to calculate a mean of R squares (default: 10)')
    optional.add_argument(
        "-f",
        "--fold",
        dest="fold",
        required=False,
        default=5,
        type=int,
        help="Specify the number of folds in a `(Stratified)KFold` (default: 5)",
    )
    optional.add_argument(
        "-n",
        "--n_permute",
        dest="n_permute",
        required=False,
        default=1000,
        type=int,
        help="The number of permutations used to calculate the p-value (default: 1000)",
    )
    optional.add_argument(
        "--predict_only",
        dest="predict_only",
        action="store_true",
        help="Only predict the risk score. Skip the permutation test.",
    )
    optional.add_argument(
        '-p',
        '--num_proc',
        dest='num_proc',
        required=False,
        type=int,
        default=1,
        help="No. worker processes for permutation (default: 1)"
    )
    optional.add_argument(
        '-S',
        '--seed',
        dest='seed',
        required=False,
        default=42,
        type=int,
        help="Seed of random state (default: 42)."
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    return result

def dawn() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas dawn",
        description="CWAS DAWN analysis",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-e",
        "--eig_vector",
        dest="eig_vector_file",
        required=True,
        type=Path,
        help="File path of eigen vectors file resulted from effective number test (*eig_vecs*.txt.gz)."
    )
    required.add_argument(
        "-c",
        "--corr_mat",
        dest="corr_mat_file",
        type=Path,
        help="File path of category correlation matrix file resulted from categorization (*correlation_matrix*.pkl)."
    )
    required.add_argument(
        "-P",
        "--permut_test",
        dest="permut_test_file",
        required=True,
        type=Path,
        help="File path of permutation test file resulted from permutation test (*permutation_test*.txt.gz)."
    )
    required.add_argument(
        "-c_count",
        "--cat_count",
        dest="category_count_file",
        required=True,
        type=Path,
        help="File path of category counts file resulted from burden test (for each variant) or sign test (for each sample).",
    )
    optional.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        default=1,
        type=int,
        help="Number of worker processes for the DAWN (default: 1).",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved (default: {}).".format(default_workspace),
    )
    optional.add_argument(
        "-r",
        "--range",
        dest="k_range",
        required=False,
        default='2,100',
        type=str,
        help="Range from start and end to find K for k-means clustering, and start must be above 1.\n"\
             "Start and end should be comma-separated. -r, and -k are mutually exclusive.\n"\
             "If you want to use range (-r) to find optimal K, -k must be None (default: 2,100).",
    )
    optional.add_argument(
        "-k",
        "--k_val",
        dest="k_val",
        required=False,
        type=int,
        default=None,
        help="K for k-means clustering. -r, and -k are mutually exclusive.\n"\
             "If `-k` is given, the value of -r is ignored (default: None).",
    )
    optional.add_argument(
        "-s",
        "--seed",
        dest="seed",
        required=False,
        default=42,
        type=int,
        help="Seed value for t-SNE (default: 42).",
    )
    optional.add_argument(
        '-T',
        '--tsen_method',
        dest='tsne_method',
        required=False,
        type=str,
        default='exact',
        choices=['barnes_hut','exact'],
        help="Gradient calculation algorithm for t-SNE, which is used in TSNE of sklearn (default: exact).\n"\
             "If the dataset is large, 'barnes_hut' is recommended."
    )
    optional.add_argument(
        "-t",
        "--tag",
        dest="tag",
        required=False,
        default=None,
        type=str,
        help="Tag used for the name of output files (e.g. intergenic, coding etc.) (default: None).",
    )
    optional.add_argument(
        "-l",
        "--lambda",
        dest="lambda_val",
        required=False,
        default=5.25,
        type=float,
        help="Lambda value for parameter tuning (default: 5.25).",
    )
    optional.add_argument(
        "-C",
        "--count_threshold",
        dest="count_threshold",
        required=False,
        type=int,
        default=20,
        help="The threshold of variant (or sample) counts, which is the least amount of variants a category should have (default: 20).",
    )
    optional.add_argument(
        "-R",
        "--corr_threshold",
        dest="corr_threshold",
        required=False,
        type=float,
        default=0.12,
        help="The threshold of correlation values between clusters. Computed by the mean value of correlation values of categories within a cluster (default: 0.12).",
    )
    optional.add_argument(
        "-S",
        "--size_threshold",
        dest="size_threshold",
        required=False,
        type=int,
        default=2,
        help="The threshold of the number of categories per cluster. The least amount of categories a cluster should have (default: 2).",
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )    
    return result


def correlation() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="cwas correlation",
        description="CWAS correlation calculation",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )
    default_workspace = dotenv.dotenv_values(dotenv_path=Path.home() / ".cwas_env").get("CWAS_WORKSPACE")
    required = result.add_argument_group("Required arguments")
    optional = result.add_argument_group("Optional arguments")
    other = result.add_argument_group("Other")
    required.add_argument(
        "-i",
        "--input_file",
        dest="cat_path",
        required=True,
        type=Path,
        help="Categorized file (*.zarr) resulted from categorization step.",
    )
    required.add_argument(
        "-v",
        "--annotated_vcf",
        dest="annot_path",
        required=False,
        type=Path,
        help="Annotated VCF file. Required for variant-level correlation matrix (--cm variant).",
    )
    optional.add_argument(
        "-o_dir",
        "--output_directory",
        dest="output_dir_path",
        required=False,
        default=default_workspace,
        type=Path,
        help="Directory where output file will be saved. (default: {})".format(default_workspace),
    )
    optional.add_argument(
        "-p",
        "--num_proc",
        dest="num_proc",
        required=False,
        type=int,
        help="Number of worker processes for the categorization. Recommend using only one-third of available cores to prevent memory errors. (default: 1)",
        default=1,
    )
    required.add_argument(
        "-cm",
        "--corr_matrix",
        dest="generate_corr_matrix",
        required=True,
        choices = ['variant','sample'],
        help="Generate a correlation matrix bewteen categories.\n * variant: use the number of variants\n * sample: use the number of samples with variants",
    )
    optional.add_argument(
        "-im",
        "--intersection_matrix",
        dest="generate_inter_matrix",
        required=False,
        action="store_true",
        help="Generate a matrix with intersected number of variants (or samples with variants) bewteen categories.",
    )
    optional.add_argument(
        '-c_info',
        '--category_info',
        dest="category_info_path",
        required=False,
        default=None,
        type=Path,
        help="Path to a text file with category information (*.category_info.txt).",
    )
    optional.add_argument(
        '-d',
        '--domain_list',
        dest="domain_list",
        required=False,
        default='all',
        type=str,
        help="Domain list to filter categories based on GENCODE domain. (default: all).\n"\
             "Available options: run_all,all,coding,noncoding,ptv,missense,damaging_missense,promoter,noncoding_wo_promoter,intron,intergenic,utr,lincRNA",
    )
    other.add_argument(
        '-h',
        '--help',
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    return result
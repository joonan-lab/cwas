import argparse
from abc import abstractmethod
from pathlib import Path
from typing import Optional
import polars as pl
import numpy as np
import pandas as pd
import re
import zarr

from cwas.core.categorization.category import Category
from cwas.core.common import cmp_two_arr
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file, check_is_dir
from cwas.utils.log import print_arg, print_progress



class BurdenTest(Runnable):
    def __init__(self, args: Optional[argparse.Namespace] = None):
        super().__init__(args)
        self._sample_info = None
        self._adj_factor = None
        self._categorization_result = None
        self._result = None
        self._phenotypes = None
        self._case_variant_cnt = None
        self._ctrl_variant_cnt = None
        self._case_carrier_cnt = None
        self._ctrl_carrier_cnt = None
        self._count_thres = None

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg("Categorized file", args.cat_path)
        print_arg("Output directory", args.output_dir_path)
        print_arg("Sample information file", args.sample_info_path)
        print_arg("Adjustment factor list", args.adj_factor_path)
        print_arg("If the number of carrier is used for burden test or not", args.use_n_carrier)

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_is_dir(args.cat_path)
        check_is_dir(args.output_dir_path)
        check_is_file(args.sample_info_path)
        if args.adj_factor_path is not None:
            check_is_file(args.adj_factor_path)

    @property
    def cat_path(self) -> Path:
        return self.args.cat_path.resolve()

    @property
    def output_dir_path(self) -> Path:
        return self.args.output_dir_path.resolve()

    @property
    def sample_info_path(self) -> Path:
        return self.args.sample_info_path.resolve()

    @property
    def adj_factor_path(self) -> Optional[Path]:
        return (
            self.args.adj_factor_path.resolve()
            if self.args.adj_factor_path
            else None
        )

    @property
    def result_path(self) -> Path:
        f_name = re.sub(r'categorization_result\.zarr\.gz|categorization_result\.zarr', 'burden_test.txt', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/"
            f"{f_name}"
        )

    @property
    def counts_path(self) -> Path:
        f_name = re.sub(r'categorization_result\.zarr\.gz|categorization_result\.zarr', 'category_counts.txt', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/"
            f"{f_name}"
        )

    @property
    def cat_info_path(self) -> Path:
        f_name = re.sub(r'categorization_result\.zarr\.gz|categorization_result\.zarr', 'category_info.txt', self.cat_path.name)
        return Path(
            f"{self.output_dir_path}/"
            f"{f_name}"
        )

    @property
    def sample_info(self) -> pd.DataFrame:
        if self._sample_info is None:
            self._sample_info = pd.read_table(
                self.sample_info_path, index_col="SAMPLE", dtype={"SAMPLE": str}, sep="\t"
            )
        return self._sample_info

    @property
    def adj_factor(self) -> pd.DataFrame:
        if self._adj_factor is None and self.adj_factor_path:
            self._adj_factor = pd.read_table(
                self.adj_factor_path, index_col="SAMPLE", dtype={"SAMPLE": str}, sep="\t"
            )
        return self._adj_factor

    @property
    def use_n_carrier(self) -> bool:
        return self.args.use_n_carrier

    @property
    def tag(self) -> str:
        return self.args.tag

    @property
    def eff_test(self) -> int:
        return self.args.eff_test

    @property
    def marker_size(self) -> float:
        return self.args.marker_size

    @property
    def font_size(self) -> float:
        return self.args.font_size

    @property
    def plot_size(self) -> float:
        return self.args.plot_size
    
    @property
    def plot_title(self) -> float:
        return self.args.plot_title

    @property
    def categorization_result(self) -> pd.DataFrame:
        if self._categorization_result is None:
            print_progress("Load the categorization result")
            root = zarr.open(self.cat_path, mode='r')
            self._categorization_result = pd.DataFrame(data=root['data'],
                              index=root['metadata'].attrs['sample_id'],
                              columns=root['metadata'].attrs['category'])
            self._categorization_result.index.name = 'SAMPLE'
            
            self.save_counts_table(form = 'raw')
            if self.adj_factor is not None:
                self._adjust_categorization_result()
        return self._categorization_result

    def _adjust_categorization_result(self):
        if not _contain_same_index(
            self._categorization_result, self.adj_factor
        ):
            raise ValueError(
                "The sample IDs from the adjustment factor list are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
        adj_factors = [
            self.adj_factor.to_dict()["AdjustFactor"][sample_id]
            for sample_id in self._categorization_result.index.values
        ]
        self._categorization_result = self._categorization_result.multiply(
            adj_factors, axis="index"
        )

    @property
    def phenotypes(self) -> np.ndarray:
        if self._phenotypes is None:
            self._phenotypes = np.array(self.sample_info.loc[self.categorization_result.index, 'PHENOTYPE'])
            
        return self._phenotypes

    @property
    def case_cnt(self) -> int:
        return (self.phenotypes == "case").sum()

    @property
    def ctrl_cnt(self) -> int:
        return (self.phenotypes == "ctrl").sum()

    @property
    def case_variant_cnt(self) -> np.ndarray:
        if self._case_variant_cnt is None:
            self._case_variant_cnt = self.categorization_result.values[
                self.phenotypes == "case", :
            ].sum(axis=0)
        return self._case_variant_cnt

    @property
    def ctrl_variant_cnt(self) -> np.ndarray:
        if self._ctrl_variant_cnt is None:
            self._ctrl_variant_cnt = self.categorization_result.values[
                self.phenotypes == "ctrl", :
            ].sum(axis=0)
        return self._ctrl_variant_cnt

    @property
    def case_carrier_cnt(self) -> np.ndarray:
        if self._case_carrier_cnt is None:
            var_counts = self.categorization_result.values
            is_carrier = np.where(var_counts > 0, 1, 0)
            self._case_carrier_cnt = is_carrier[
                self.phenotypes == "case", :
            ].sum(axis=0)
        return self._case_carrier_cnt

    @property
    def ctrl_carrier_cnt(self) -> np.ndarray:
        if self._ctrl_carrier_cnt is None:
            var_counts = self.categorization_result.values
            is_carrier = np.where(var_counts > 0, 1, 0)
            self._ctrl_carrier_cnt = is_carrier[
                self.phenotypes == "ctrl", :
            ].sum(axis=0)
        return self._ctrl_carrier_cnt

    @property
    def category_table(self) -> pd.DataFrame:
        categories = [
            Category.from_str(category_str).to_dict()
            for category_str in self.categorization_result.columns.values
        ]
        return pd.DataFrame(
            categories, index=self.categorization_result.columns.values
        )

    def run(self):
        if not _contain_same_index(
            self.categorization_result, self.sample_info
        ):
            raise ValueError(
                "The sample IDs from the sample information are "
                "not the same with the sample IDs "
                "from the categorization result."
            )
        if self.use_n_carrier:
            self.count_carrier_for_each_category()
            self.calculate_relative_risk_with_n_carrier()
        else:
            self.count_variant_for_each_category()
            self.calculate_relative_risk()        
        self.run_burden_test()
        self.concat_category_info()
        self.save_result()
        self.save_counts_table(form = 'adj')
        self.save_category_info()
        self.update_env()

    def count_variant_for_each_category(self):
        print_progress("Count the number of variants for each category")
        variant_cnt_arr = np.concatenate(
            [
                self.case_variant_cnt[:, np.newaxis],
                self.ctrl_variant_cnt[:, np.newaxis],
            ],
            axis=1,
        )
        self._result = pd.DataFrame(
            variant_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_DNV_Count", "Ctrl_DNV_Count"],
        )

    def count_carrier_for_each_category(self):
        print_progress("Count the number of carriers in each category")
        carrier_cnt_arr = np.concatenate(
            [
                self.case_carrier_cnt[:, np.newaxis],
                self.ctrl_carrier_cnt[:, np.newaxis],
            ],
            axis=1,
        )
        self._result = pd.DataFrame(
            carrier_cnt_arr,
            index=self.categorization_result.columns.values,
            columns=["Case_Carrier_Count", "Ctrl_Carrier_Count"],
        )

    def calculate_relative_risk(self):
        print_progress("Calculate relative risks for each category")
        np.seterr(divide='ignore')
        normalized_case_variant_cnt = self.case_variant_cnt / self.case_cnt
        normalized_ctrl_variant_cnt = self.ctrl_variant_cnt / self.ctrl_cnt
        self._result["Relative_Risk"] = (
            normalized_case_variant_cnt / normalized_ctrl_variant_cnt
        )

    def calculate_relative_risk_with_n_carrier(self):
        print_progress("Calculate relative risks for each category")
        case_carrier_rate = self.case_carrier_cnt / self.case_cnt
        ctrl_carrier_rate = self.ctrl_carrier_cnt / self.ctrl_cnt
        self._result["Relative_Risk"] = (
            case_carrier_rate / ctrl_carrier_rate
        )

    @abstractmethod
    def run_burden_test(self):
        raise RuntimeError(
            "This method cannot be called via the instance of BurdenTest."
        )

    def concat_category_info(self):
        self._result = pd.concat([self.category_table, self._result], axis=1)
        self._result.index.name = "Category"

    def save_result(self):
        print_progress(f"Save the result to the file {self.result_path}")
        self._result.to_csv(self.result_path, sep="\t")

    def update_env(self):
        self.set_env("BURDEN_TEST_RESULT", self.result_path)
        self.save_env()
        
    def save_counts_table(self, form: str):
        if form == 'raw':
            if self.use_n_carrier:
                self._raw_counts = pd.DataFrame({'Raw_counts': (self.categorization_result > 0).sum(axis=0)},
                                                index=self.categorization_result.columns)
                self._raw_counts.index.name = 'Category'

            else:
                self._raw_counts = pd.DataFrame({'Raw_counts': self.categorization_result.sum(axis=0)},
                                                index= self.categorization_result.sum(axis=0).index)
                self._raw_counts.index.name = 'Category'
        elif form =='adj':
            if self.use_n_carrier:
                self._adj_counts = pd.DataFrame({'Adj_counts': self._result['Case_Carrier_Count'] + self._result['Ctrl_Carrier_Count']})
            else:
                self._adj_counts = pd.DataFrame({'Adj_counts': self._result['Case_DNV_Count'] + self._result['Ctrl_DNV_Count']})
            self._counts_table = pd.merge(self._raw_counts, self._adj_counts, on='Category')
            print_progress(f"Save the category counts to the file {self.counts_path}")
            self._counts_table.to_csv(self.counts_path, sep="\t")
    
    def save_category_info(self):
        cat_set = self._result.loc[:, ['variant_type', 'gene_set', 'functional_score', 'gencode', 'functional_annotation']]
        cat_set = apply_region_mapping(cat_set)
        print_progress(f"Save the category information to the file {self.cat_info_path}")
        cat_set.to_csv(self.cat_info_path, sep="\t")



def _contain_same_index(table1: pd.DataFrame, table2: pd.DataFrame) -> bool:
    return cmp_two_arr(table1.index.values, table2.index.values)

def apply_region_mapping(df: pd.DataFrame):
    coding_region = ['CodingRegion', 'MissenseRegion', 'SilentRegion', 'PTVRegion',
                     'InFrameRegion', 'DamagingMissenseRegion', 'FrameshiftRegion']
    noncoding_region = ['NoncodingRegion', 'IntronRegion', 'lincRnaRegion',
                        'IntergenicRegion', 'OtherTranscriptRegion', 'PromoterRegion',
                        'UTRsRegion', 'SpliceSiteNoncanonRegion']
    
    region_mapping = {
        'is_coding': lambda x: x['gencode'].isin(coding_region).astype(int),
        'is_coding_no_ptv': lambda x: x['gencode'].isin(set(coding_region) - set(['CodingRegion', 'PTVRegion', 'FrameshiftRegion'])).astype(int),
        'is_PTV': lambda x: (x['gencode'] == 'PTVRegion').astype(int),
        'is_missense': lambda x: (x['gencode'] == 'MissenseRegion').astype(int),
        'is_damaging_missense': lambda x: (x['gencode'] == 'DamagingMissenseRegion').astype(int),
        'is_noncoding': lambda x: x['gencode'].isin(noncoding_region).astype(int),
        'is_noncoding_wo_promoter': lambda x: x['gencode'].isin(set(noncoding_region) - set(['NoncodingRegion', 'PromoterRegion'])).astype(int),
        'is_promoter': lambda x: (x['gencode'] == 'PromoterRegion').astype(int),
        'is_intron': lambda x: (x['gencode'] == 'IntronRegion').astype(int),
        'is_intergenic': lambda x: (x['gencode'] == 'IntergenicRegion').astype(int),
        'is_UTR': lambda x: (x['gencode'] == 'UTRsRegion').astype(int),
        'is_lincRNA': lambda x: ((x['gene_set'] == 'lincRNA') | (x['gencode'] == 'lincRnaRegion')).astype(int)
    }

    for col, condition in region_mapping.items():
        df[col] = condition(df)
    
    return df
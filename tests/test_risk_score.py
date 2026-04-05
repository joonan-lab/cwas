"""
Tests for cwas.risk_score module
"""
import pytest
import argparse
import numpy as np
from pathlib import Path
from unittest.mock import patch
from collections import defaultdict
from cwas.risk_score import RiskScore


@pytest.fixture
def make_risk_score():
    def _make(**overrides):
        defaults = dict(
            sample_info_path=Path("/tmp/sample_info.txt"),
            categorization_result_path=Path("/tmp/categorization.zarr"),
            adj_factor_path=None,
            category_set_path=Path("/tmp/category_set.txt"),
            domain_list="all",
            tag="test",
            use_n_carrier=False,
            do_each_one=False,
            leave_one_out=False,
            ctrl_thres=5,
            train_set_f=0.8,
            num_reg=3,
            fold=5,
            n_permute=10,
            predict_only=False,
            num_proc=1,
            seed=42,
            plotsize="8,6",
            fontsize=12.0,
            feature_selection_group="gene_set",
            output_dir_path=Path("/tmp/output"),
        )
        defaults.update(overrides)
        args = argparse.Namespace(**defaults)
        with patch.object(RiskScore, '_print_args'), \
             patch.object(RiskScore, '_check_args_validity'), \
             patch('cwas.risk_score.importr'):
            return RiskScore(args)
    return _make


class TestRiskScoreInit:
    def test_init_creates_instance(self, make_risk_score):
        rs = make_risk_score()
        assert rs._sample_info is None
        assert rs._categorization_result is None
        assert isinstance(rs._result_dict, defaultdict)

    def test_tag_property(self, make_risk_score):
        rs = make_risk_score(tag="mytest")
        assert rs.tag == "mytest"

    def test_use_n_carrier_property(self, make_risk_score):
        rs = make_risk_score(use_n_carrier=True)
        assert rs.use_n_carrier is True

    def test_ctrl_thres_property(self, make_risk_score):
        rs = make_risk_score(ctrl_thres=10)
        assert rs.ctrl_thres == 10

    def test_train_set_f_property(self, make_risk_score):
        rs = make_risk_score(train_set_f=0.7)
        assert rs.train_set_f == 0.7

    def test_num_reg_property(self, make_risk_score):
        rs = make_risk_score(num_reg=5)
        assert rs.num_reg == 5

    def test_fold_property(self, make_risk_score):
        rs = make_risk_score(fold=10)
        assert rs.fold == 10

    def test_seed_property(self, make_risk_score):
        rs = make_risk_score(seed=123)
        assert rs.seed == 123

    def test_num_proc_property(self, make_risk_score):
        rs = make_risk_score(num_proc=4)
        assert rs.num_proc == 4

    def test_domain_list_all(self, make_risk_score):
        rs = make_risk_score(domain_list="all")
        assert rs.domain_list == ['all']


class TestFeatureSelectionGroup:
    def test_single_value(self, make_risk_score):
        rs = make_risk_score(feature_selection_group="gene_set")
        assert rs.feature_selection_group == ["gene_set"]

    def test_multiple_values(self, make_risk_score):
        rs = make_risk_score(feature_selection_group="gene_set, functional_score")
        result = rs.feature_selection_group
        assert "gene_set" in result
        assert "functional_score" in result

    def test_invalid_value(self, make_risk_score):
        rs = make_risk_score(feature_selection_group="invalid_feature")
        with pytest.raises(ValueError, match="Invalid feature selection group"):
            _ = rs.feature_selection_group

    def test_case_insensitive(self, make_risk_score):
        rs = make_risk_score(feature_selection_group="GENE_SET")
        assert rs.feature_selection_group == ["gene_set"]


class TestCustomCVFolds:
    def test_basic(self, make_risk_score):
        rs = make_risk_score(fold=3)
        foldid = rs._custom_cv_folds(nobs=30, seed=42)
        assert len(foldid) == 30
        assert np.max(foldid) <= 3
        assert np.min(foldid) >= 0

    def test_deterministic(self, make_risk_score):
        rs = make_risk_score(fold=5)
        fold1 = rs._custom_cv_folds(nobs=50, seed=42)
        fold2 = rs._custom_cv_folds(nobs=50, seed=42)
        np.testing.assert_array_equal(fold1, fold2)

    def test_different_seeds(self, make_risk_score):
        rs = make_risk_score(fold=5)
        fold1 = rs._custom_cv_folds(nobs=50, seed=42)
        fold2 = rs._custom_cv_folds(nobs=50, seed=123)
        assert not np.array_equal(fold1, fold2)


class TestCheckDomainList:
    def test_valid(self, make_risk_score):
        rs = make_risk_score()
        result = rs._check_domain_list("missense", ["Missense", "LOF"])
        assert result == "Missense"

    def test_invalid(self, make_risk_score):
        rs = make_risk_score()
        with pytest.raises(ValueError, match="Invalid domain name"):
            rs._check_domain_list("invalid", ["Missense", "LOF"])

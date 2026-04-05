"""
Tests for cwas.core.configuration.settings module
"""
import pytest
from cwas.core.configuration.settings import (
    get_default_domains,
    get_domain_types,
    get_redundant_domain_pairs,
)


class TestGetDefaultDomains:
    def test_returns_dict(self):
        result = get_default_domains()
        assert isinstance(result, dict)

    def test_contains_required_keys(self):
        result = get_default_domains()
        required = ["variant_type", "functional_score", "gene_set", "gencode", "functional_annotation"]
        for key in required:
            assert key in result

    def test_variant_type_values(self):
        result = get_default_domains()
        assert result["variant_type"] == ["All", "SNV", "Indel"]

    def test_functional_score_contains_all(self):
        result = get_default_domains()
        assert "All" in result["functional_score"]

    def test_gene_set_contains_any(self):
        result = get_default_domains()
        assert "Any" in result["gene_set"]

    def test_gencode_starts_with_any(self):
        result = get_default_domains()
        assert result["gencode"][0] == "Any"

    def test_gencode_has_coding_regions(self):
        result = get_default_domains()
        expected = ["CodingRegion", "PTVRegion", "MissenseRegion", "NoncodingRegion", "IntronRegion"]
        for val in expected:
            assert val in result["gencode"]

    def test_returns_independent_copies(self):
        result1 = get_default_domains()
        result1["variant_type"].append("Modified")
        result2 = get_default_domains()
        assert "Modified" not in result2["variant_type"]

    def test_all_values_are_lists(self):
        result = get_default_domains()
        for key, value in result.items():
            assert isinstance(value, list)


class TestGetDomainTypes:
    def test_returns_list(self):
        result = get_domain_types()
        assert isinstance(result, list)

    def test_has_five_types(self):
        result = get_domain_types()
        assert len(result) == 5

    def test_contains_expected_types(self):
        result = get_domain_types()
        expected = ["variant_type", "functional_score", "gene_set", "gencode", "functional_annotation"]
        for t in expected:
            assert t in result

    def test_returns_independent_copies(self):
        result1 = get_domain_types()
        result1.append("NewDomain")
        result2 = get_domain_types()
        assert "NewDomain" not in result2


class TestGetRedundantDomainPairs:
    def test_returns_dict(self):
        result = get_redundant_domain_pairs()
        assert isinstance(result, dict)

    def test_has_two_pair_types(self):
        result = get_redundant_domain_pairs()
        assert len(result) == 2

    def test_variant_type_gencode_pair_exists(self):
        result = get_redundant_domain_pairs()
        assert ("variant_type", "gencode") in result

    def test_gene_set_gencode_pair_exists(self):
        result = get_redundant_domain_pairs()
        assert ("gene_set", "gencode") in result

    def test_variant_type_gencode_has_expected_pairs(self):
        result = get_redundant_domain_pairs()
        pairs = result[("variant_type", "gencode")]
        assert ("All", "FrameshiftRegion") in pairs
        assert ("SNV", "FrameshiftRegion") in pairs
        assert ("Indel", "DamagingMissenseRegion") in pairs

    def test_gene_set_gencode_has_expected_pairs(self):
        result = get_redundant_domain_pairs()
        pairs = result[("gene_set", "gencode")]
        assert ("Any", "CodingRegion") in pairs
        assert ("Any", "lincRnaRegion") in pairs

    def test_returns_independent_copies(self):
        result1 = get_redundant_domain_pairs()
        key = ("variant_type", "gencode")
        result1[key].add(("New", "Pair"))
        result2 = get_redundant_domain_pairs()
        assert ("New", "Pair") not in result2[key]

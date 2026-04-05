"""
Tests for cwas.core.categorization.categorizer module
"""
import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch
from cwas.core.categorization.categorizer import Categorizer


class TestCategorizerInit:
    def test_init_basic(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any", "CodingRegion"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        gene_matrix = {"BRCA1": {"Any"}, "BRCA2": {"Any"}}
        categorizer = Categorizer(category_domain, gene_matrix, "Score", 0.5)
        assert categorizer._category_domain == category_domain
        assert categorizer._gene_matrix == gene_matrix
        assert categorizer._mis_info_key == "Score"
        assert categorizer._mis_thres == 0.5

    def test_init_empty_gene_matrix(self):
        category_domain = {"variant_type": ["All"], "gene_set": ["Any"]}
        categorizer = Categorizer(category_domain, {}, "LoF", 0.3)
        assert categorizer._gene_matrix == {}


class TestAnnotateVariantType:
    def test_snv_detection(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        vcf_data = pd.DataFrame({"REF": ["A"], "ALT": ["T"]})
        result = categorizer.annotate_variant_type(vcf_data)
        # All(0) + SNV(1) = 2^0 + 2^1 = 3
        assert result[0] == 3

    def test_indel_insertion(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        vcf_data = pd.DataFrame({"REF": ["A"], "ALT": ["ATT"]})
        result = categorizer.annotate_variant_type(vcf_data)
        # All(0) + Indel(2) = 2^0 + 2^2 = 5
        assert result[0] == 5

    def test_indel_deletion(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        vcf_data = pd.DataFrame({"REF": ["ATT"], "ALT": ["A"]})
        result = categorizer.annotate_variant_type(vcf_data)
        assert result[0] == 5

    def test_multiple_variants(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        vcf_data = pd.DataFrame({
            "REF": ["A", "A", "ATT"],
            "ALT": ["T", "ATT", "A"]
        })
        result = categorizer.annotate_variant_type(vcf_data)
        assert len(result) == 3
        assert result[0] == 3  # SNV
        assert result[1] == 5  # Indel
        assert result[2] == 5  # Indel


class TestParseAnnotationInt:
    def test_single_bit(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        result = categorizer.parse_annotation_int(2, "variant_type")
        assert "SNV" in result

    def test_multiple_bits(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        result = categorizer.parse_annotation_int(3, "variant_type")
        assert "All" in result
        assert "SNV" in result

    def test_zero(self):
        category_domain = {
            "variant_type": ["All", "SNV", "Indel"],
            "gene_set": ["Any"],
            "gencode": ["Any"],
            "functional_score": ["All"],
            "functional_annotation": ["Any"]
        }
        categorizer = Categorizer(category_domain, {}, "Score", 0.5)
        result = categorizer.parse_annotation_int(0, "variant_type")
        assert len(result) == 0

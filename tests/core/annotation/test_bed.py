"""
Tests for cwas.core.annotation.bed module
"""
import pytest
from unittest.mock import MagicMock, patch
from cwas.core.annotation.bed import annotate


class TestAnnotateInit:
    def test_init_basic(self):
        ann = annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=4)
        assert ann.in_vcf_gz_path == "/in.vcf.gz"
        assert ann.out_vcf_path == "/out.vcf"
        assert ann.annot_bed_path == "/annot.bed.gz"
        assert ann.num_proc == 4

    def test_num_proc_capped_at_22(self):
        ann = annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=50)
        assert ann.num_proc == 22

    def test_num_proc_not_capped_below_23(self):
        ann = annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=10)
        assert ann.num_proc == 10

    def test_num_proc_boundary_22(self):
        ann = annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=22)
        assert ann.num_proc == 22

    def test_num_proc_boundary_23(self):
        ann = annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=23)
        assert ann.num_proc == 22


class TestChrAnnotate:
    @pytest.fixture
    def annotator(self):
        return annotate("/in.vcf.gz", "/out.vcf", "/annot.bed.gz", num_proc=1)

    def test_snv_no_bed_overlap(self, annotator):
        mock_vcf = MagicMock()
        mock_bed = MagicMock()
        variant = ("chr1", "100", ".", "A", "T")
        mock_vcf.fetch.return_value = iter([variant, None])
        mock_bed.fetch.return_value = iter([])
        with patch('cwas.core.annotation.bed.pysam.TabixFile') as mock_tabix:
            mock_tabix.side_effect = [mock_vcf, mock_bed]
            result = annotator.chr_annotate("chr1")
        assert len(result) == 1
        assert "ANNOT=0" in result[0]

    def test_snv_with_bed_overlap(self, annotator):
        mock_vcf = MagicMock()
        mock_bed = MagicMock()
        variant = ("chr1", "100", ".", "A", "T")
        bed_region = ("chr1", "90", "110", "5")
        mock_vcf.fetch.return_value = iter([variant, None])
        mock_bed.fetch.return_value = iter([bed_region])
        with patch('cwas.core.annotation.bed.pysam.TabixFile') as mock_tabix:
            mock_tabix.side_effect = [mock_vcf, mock_bed]
            result = annotator.chr_annotate("chr1")
        assert len(result) == 1
        assert "ANNOT=5" in result[0]

    def test_bitwise_annotation_combination(self, annotator):
        mock_vcf = MagicMock()
        mock_bed = MagicMock()
        variant = ("chr1", "100", ".", "A", "T")
        bed1 = ("chr1", "90", "110", "1")
        bed2 = ("chr1", "95", "115", "2")
        mock_vcf.fetch.return_value = iter([variant, None])
        mock_bed.fetch.return_value = iter([bed1, bed2])
        with patch('cwas.core.annotation.bed.pysam.TabixFile') as mock_tabix:
            mock_tabix.side_effect = [mock_vcf, mock_bed]
            result = annotator.chr_annotate("chr1")
        assert len(result) == 1
        assert "ANNOT=3" in result[0]  # 1 | 2 = 3

    def test_multiple_variants(self, annotator):
        mock_vcf = MagicMock()
        mock_bed = MagicMock()
        variants = [
            ("chr1", "100", ".", "A", "T"),
            ("chr1", "200", ".", "G", "C"),
        ]
        mock_vcf.fetch.return_value = iter(variants + [None])
        mock_bed.fetch.return_value = iter([])
        with patch('cwas.core.annotation.bed.pysam.TabixFile') as mock_tabix:
            mock_tabix.side_effect = [mock_vcf, mock_bed]
            result = annotator.chr_annotate("chr1")
        assert len(result) == 2

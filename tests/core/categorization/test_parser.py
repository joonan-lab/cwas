from io import StringIO

from cwas.core.categorization import parser
from pandas import Series


def test_parse_info_field():
    csq_field_names = ["Test1", "Test2", "Test3"]
    info_field = (
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="'
        f"Consequence annotations from Ensembl VEP. Format: "
        f"{'|'.join(csq_field_names)}\">"
    )
    assert parser._parse_vcf_info_field(info_field) == csq_field_names


def test_parse_annot_field():
    annot_field_names = ["ANNOT1", "ANNOT2", "ANNOT3"]
    annot_field = f"##INFO=<ID=ANNOT,Key={'|'.join(annot_field_names)}>"
    assert parser._parse_annot_field(annot_field) == annot_field_names


def test_parse_vcf_header_line():
    vcf_column_titles = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
    ]

    header_line = "\t".join(vcf_column_titles)
    assert parser._parse_vcf_header_line(f"#{header_line}") == vcf_column_titles


def test_parse_info_column():
    info_column = Series(
        [
            "SAMPLE=A;BATCH=1;CSQ=1|2|3;ANNOT=7",
            "SAMPLE=B;BATCH=1;CSQ=2|3|1;ANNOT=5",
            "SAMPLE=C;BATCH=1;CSQ=3|2|1;ANNOT=3",
            "SAMPLE=D;BATCH=1;CSQ=3|1|2;ANNOT=1",
        ]
    )
    csq_field_names = ["CSQ1", "CSQ2", "CSQ3"]
    annot_field_names = ["ANNOT1", "ANNOT2", "ANNOT3"]
    info_df = parser._parse_info_column(
        info_column, csq_field_names, annot_field_names
    )
    assert info_df["CSQ1"].to_list() == ["1", "2", "3", "3"]
    assert info_df["CSQ2"].to_list() == ["2", "3", "2", "1"]
    assert info_df["CSQ3"].to_list() == ["3", "1", "1", "2"]
    assert info_df["ANNOT1"].to_list() == [1, 1, 1, 1]
    assert info_df["ANNOT2"].to_list() == [1, 0, 1, 0]
    assert info_df["ANNOT3"].to_list() == [1, 1, 0, 0]


def test_parse_gene_matrix():
    gene_matrix_lines = [
        "\t".join(["gene_id", "gene_name", "asd1", "asd2", "asd3", "asd4"]),
        "\t".join(["id1", "gene1", "1", "0", "1", "1"]),
        "\t".join(["id2", "gene2", "0", "0", "1", "1"]),
        "\t".join(["id3", "gene3", "0", "1", "1", "0"]),
    ]
    result = parser._parse_gene_matrix(StringIO("\n".join(gene_matrix_lines)))
    assert result["gene1"] == {"asd1", "asd3", "asd4"}
    assert result["gene2"] == {"asd3", "asd4"}
    assert result["gene3"] == {"asd2", "asd3"}

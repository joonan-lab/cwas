from cwas.core.categorization import parser


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
    header_line = f"#'\t'.join(vcf_column_titles)"
    assert parser._parse_vcf_header_line(header_line) == vcf_column_titles


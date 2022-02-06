"""
Module for the categorization of the variants into CWAS categories and 
counting the number of variants in each category.
CWAS category is a combination of annotation terms 
from each of annotation groups.

There are currently 5 groups of annotation terms.

--- The groups of the annotation terms for CWAS ---
    1. Variant types (variant_type)
    2. Conservation (conservation)
    3. Gene lists (gene_list)
    4. GENCODE annotation categories (gencode)
    5. Functional annotation categories (region)

"""
import numpy as np
import pandas as pd

from itertools import product
from collections import defaultdict
from cwas.core.categorization.category import Category

# TODO: Improve the readability of this code
class Categorizer:
    def __init__(self, category_domain: dict, gene_matrix: dict) -> None:
        self._category_domain = category_domain
        self._gene_matrix = gene_matrix

    def categorize_variant(self, annotated_vcf: pd.DataFrame):
        """ Categorize the variants in the input data frame 
        into the combinations of the annotation terms in the 5 groups
        and return the dictionary that contains the distribution of the variants 
        for the combinations. These combinations are defined as CWAS categories.
        """
        # Categorize by the annotation terms for each variant
        result = defaultdict(int)

        for annotation_term_lists in self.annotate_each_variant(annotated_vcf):
            for combination in product(*annotation_term_lists):
                result[Category(*combination)] += 1

        return result

    def annotate_each_variant(self, annotated_vcf):
        """Newly annotated each variant using CWAS annotation terms"""
        # For annotating each variant with multiple annotation terms
        # from each group efficiently, "Annotation integer" is used.
        # Annotation integer: A bitwise representation of the annotation
        # of each variant where each bit means each annotation term
        # e.g. If a list of annotation terms is ['A', 'B', 'C'] and
        # the annotation integer is 0b101, it means that
        # the variant is annotated as 'A' and 'B'.
        variant_type_annot_ints = self.annotate_variant_type(annotated_vcf)
        conservation_annot_ints = annot_conservation(
            annotated_vcf, get_idx_dict(self._category_domains["conservation"])
        )
        gene_list_annot_ints = annot_gene_list(
            annotated_vcf,
            get_idx_dict(self._category_domains["gene_list"]),
            self._gene_matrix,
        )
        gencode_annot_ints = annot_gencode(
            annotated_vcf,
            get_idx_dict(self._category_domains["gencode"]),
            self._gene_matrix,
        )
        region_annot_ints = annot_region(
            annotated_vcf, get_idx_dict(self._category_domains["region"])
        )

        for (
            variant_type_annot_int,
            conservation_annot_int,
            gene_list_annot_int,
            gencode_annot_int,
            region_annot_int,
        ) in zip(
            variant_type_annot_ints,
            conservation_annot_ints,
            gene_list_annot_ints,
            gencode_annot_ints,
            region_annot_ints,
        ):
            yield (
                parse_annot_int(
                    variant_type_annot_int,
                    self._category_domains["variant_type"],
                ),
                parse_annot_int(
                    conservation_annot_int,
                    self._category_domains["conservation"],
                ),
                parse_annot_int(
                    gene_list_annot_int, self._category_domains["gene_list"]
                ),
                parse_annot_int(
                    gencode_annot_int, self._category_domains["gencode"]
                ),
                parse_annot_int(
                    region_annot_int, self._category_domains["region"]
                ),
            )

    def annotate_variant_type(self, annotated_vcf: pd.DataFrame) -> list:
        variant_type_annotation_idx = get_idx_dict(
            self._category_domains["variant_type"]
        )
        refs = annotated_vcf["REF"].values
        alts = annotated_vcf["ALT"].values

        is_snv_arr = (
            (np.vectorize(len)(refs) == 1) & (np.vectorize(len)(alts) == 1)
        ).astype(np.int32)
        annotation_int_conv = (
            lambda is_snv: 2 ** variant_type_annotation_idx["SNV"]
            if is_snv
            else 2 ** variant_type_annotation_idx["Indel"]
        )
        annotation_ints = np.vectorize(annotation_int_conv)(is_snv_arr)
        annotation_ints += 2 ** variant_type_annotation_idx["All"]

        return annotation_ints


def get_idx_dict(list_: list) -> dict:
    return {item: i for i, item in enumerate(list_)}


# Step 2: Annotate by conservationervation scores
def annot_conservation(
    variant_df: pd.DataFrame, conservation_annot_idx_dict: dict
):
    phylop_conv_func = (
        lambda x: -2.0 if x == "" else max(map(float, x.split("&")))
    )
    phast_conv_func = (
        lambda x: 0.0 if x == "" else max(map(float, x.split("&")))
    )

    phylop_scores = np.vectorize(phylop_conv_func)(
        variant_df["phyloP46wayVt"].values
    )
    phast_scores = np.vectorize(phast_conv_func)(
        variant_df["phastCons46wayVt"].values
    )

    is_phylop_conservation_arr = (phylop_scores >= 2.0).astype(np.int32)
    is_phast_conservation_arr = (phast_scores >= 0.2).astype(np.int32)

    annot_ints = np.vectorize(
        lambda is_phylop_conservation: 2
        ** conservation_annot_idx_dict["phyloP46wayVt"]
        if is_phylop_conservation
        else 0
    )(is_phylop_conservation_arr)
    annot_ints += np.vectorize(
        lambda is_phast_conservation: 2
        ** conservation_annot_idx_dict["phastCons46wayVt"]
        if is_phast_conservation
        else 0
    )(is_phast_conservation_arr)
    annot_ints += 2 ** conservation_annot_idx_dict["All"]

    return annot_ints


# Step 3: Annotate by the names of the gene lists
#         where the variant-associated genes are involved
def annot_gene_list(
    variant_df: pd.DataFrame,
    gene_annot_idx_dict: dict,
    gene_list_set_dict: dict,
):
    gene_symbols = variant_df["SYMBOL"].values
    gene_nearests = variant_df["NEAREST"].values
    gencodes = variant_df["Consequence"].values  # GENCODE annotations

    annot_int_list = []
    annot_int_dict = (
        {}
    )  # Key: a gene symbol, Value: its category integer (For memorization)

    for symbol, nearest, gencode in zip(gene_symbols, gene_nearests, gencodes):
        gene = (
            nearest
            if "downstream_gene_variant" in gencode
            or "intergenic_variant" in gencode
            else symbol
        )
        annot_int = annot_int_dict.get(gene, 0)

        if annot_int == 0:
            gene_list_set = gene_list_set_dict.get(gene, set())

            if gene_list_set:
                for gene_cat in gene_annot_idx_dict:
                    if gene_cat in gene_list_set:
                        annot_int += 2 ** gene_annot_idx_dict[gene_cat]

        annot_int_list.append(annot_int)

    annot_ints = np.asarray(annot_int_list)
    annot_ints += 2 ** gene_annot_idx_dict["Any"]

    return annot_ints


# Step 4: Annotate by the gencodes (GENCODE annotations)
def annot_gencode(
    variant_df: pd.DataFrame,
    gencode_annot_idx_dict: dict,
    gene_list_set_dict: dict,
):
    gene_symbols = variant_df["SYMBOL"].values
    gene_nearests = variant_df["NEAREST"].values
    gencodes = variant_df["Consequence"].values
    polyphens = variant_df["PolyPhen"].values

    annot_int_list = []

    for symbol, nearest, gencode, polyphen in zip(
        gene_symbols, gene_nearests, gencodes, polyphens
    ):
        gene = (
            nearest
            if "downstream_gene_variant" in gencode
            or "intergenic_variant" in gencode
            else symbol
        )
        gene_cat_set = gene_list_set_dict.get(gene, set())
        annot_int = 0
        is_in_coding = False

        if "geneSet_Protein_Coding" in gene_cat_set:
            is_in_coding = True
            annot_int += 2 ** gencode_annot_idx_dict["CodingRegion"]

            # Coding region
            if (
                "stop_gained" in gencode
                or "splice_donor" in gencode
                or "splice_acceptor" in gencode
            ):
                annot_int += 2 ** gencode_annot_idx_dict["LoFRegion"]
            elif (
                "frameshift_variant" in gencode
                or "transcript_amplification" in gencode
                or "transcript_ablation" in gencode
            ):
                annot_int += 2 ** gencode_annot_idx_dict["LoFRegion"]
                annot_int += 2 ** gencode_annot_idx_dict["FrameshiftRegion"]
            elif (
                "missense_variant" in gencode
                or "start_lost" in gencode
                or "stop_lost" in gencode
            ):
                annot_int += 2 ** gencode_annot_idx_dict["MissenseRegion"]

                if "probably_damaging" in polyphen:
                    annot_int += (
                        2 ** gencode_annot_idx_dict["MissenseHVARDRegionSimple"]
                    )

            elif (
                "inframe_deletion" in gencode or "inframe_insertion" in gencode
            ):
                annot_int += 2 ** gencode_annot_idx_dict["InFrameRegion"]
            elif "synonymous_variant" in gencode:
                annot_int += 2 ** gencode_annot_idx_dict["SilentRegion"]
            elif (
                "stop_retained_variant" not in gencode
                and "incomplete_terminal_codon_variant" not in gencode
                and "protein_altering_variant" not in gencode
                and "coding_sequence_variant" not in gencode
            ):
                # Noncoding
                annot_int = 0
                is_in_coding = False

        if not is_in_coding:
            annot_int += 2 ** gencode_annot_idx_dict["NoncodingRegion"]

            if "_UTR_" in gencode:
                annot_int += 2 ** gencode_annot_idx_dict["UTRsRegion"]
            elif "upstream_gene_variant" in gencode:
                annot_int += 2 ** gencode_annot_idx_dict["PromoterRegion"]
            elif "intron_variant" in gencode:
                annot_int += 2 ** gencode_annot_idx_dict["IntronRegion"]
            elif "splice_region_variant" in gencode:
                annot_int += (
                    2 ** gencode_annot_idx_dict["SpliceSiteNoncanonRegion"]
                )
            elif (
                "downstream_gene_variant" in gencode
                or "intergenic_variant" in gencode
            ):
                annot_int += 2 ** gencode_annot_idx_dict["IntergenicRegion"]
            elif "geneSet_Protein_Coding" not in gene_cat_set:
                if "geneSet_Antisense" in gene_cat_set:
                    annot_int += 2 ** gencode_annot_idx_dict["AntisenseRegion"]
                elif "geneSet_lincRNA" in gene_cat_set:
                    annot_int += 2 ** gencode_annot_idx_dict["lincRnaRegion"]
                else:
                    annot_int += (
                        2 ** gencode_annot_idx_dict["OtherTranscriptRegion"]
                    )

        annot_int_list.append(annot_int)

    annot_ints = np.asarray(annot_int_list)
    annot_ints += 2 ** gencode_annot_idx_dict["Any"]

    return annot_ints


# Step 5: Annotate by the annotation terms of the regions where the variants are
def annot_region(variant_df: pd.DataFrame, region_annot_idx_dict: dict):
    annot_ints = np.zeros(len(variant_df.index), dtype=int)

    for region in region_annot_idx_dict:
        if region == "Any":
            continue

        region_vals = variant_df[region].values.astype(np.int32)
        annot_int_conv_func = lambda x: 2 ** region_annot_idx_dict[region] * x
        annot_ints += np.vectorize(annot_int_conv_func)(region_vals)

    annot_ints += 2 ** region_annot_idx_dict["Any"]

    return annot_ints


def parse_annot_int(annot_int: int, all_annot_terms: list) -> list:
    """ From the list of annotation terms, 
    choose the appropriate subset by parsing the annotation integer.
    """
    annot_idx = 0
    annot_terms = []

    while annot_int != 0:
        if annot_int % 2 == 1:
            annot_terms.append(all_annot_terms[annot_idx])

        annot_int >>= 1
        annot_idx += 1

    return annot_terms

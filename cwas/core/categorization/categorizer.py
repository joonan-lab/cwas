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
from collections import defaultdict
from itertools import product

import numpy as np
import pandas as pd
from cwas.core.categorization.category import Category
from cwas.core.categorization.utils import extract_sublist_by_int, get_idx_dict


class Categorizer:
    def __init__(self, category_domain: dict, gene_matrix: dict) -> None:
        self._category_domain = category_domain
        self._gene_matrix = gene_matrix
        self._default_phylop = -2.0
        self._default_phastcons = 0.0
        self._phylop_cutoff = 2.0
        self._phastcons_cutoff = 0.2
        self._default_jarvis099 = 0.0
        self._default_constraintZ4 = 0.0
        self._jarvis099_cutoff = 0.99
        self._constraintZ4_cutoff = 4.0
        self._default_Splice02max = 0.0
        self._default_Splice05max = 0.0
        self._default_Splice08max = 0.0
        self._Splice02max_cutoff = 0.2
        self._Splice05max_cutoff = 0.5
        self._Splice08max_cutoff = 0.8

    @property
    def phylop_cutoff(self):
        return self._phylop_cutoff

    @phylop_cutoff.setter
    def phylop_cutoff(self, value):
        self._phylop_cutoff = value

    @property
    def phastcons_cutoff(self):
        return self._phastcons_cutoff

    @phastcons_cutoff.setter
    def phastcons_cutoff(self, value):
        self._phastcons_cutoff = value

    @property
    def jarvis099_cutoff(self):
        return self._jarvis099_cutoff

    @jarvis099_cutoff.setter
    def jarvis099_cutoff(self, value):
        self._jarvis099_cutoff = value

    @property
    def constraintZ4_cutoff(self):
        return self._constraintZ4_cutoff

    @constraintZ4_cutoff.setter
    def constraintZ4_cutoff(self, value):
        self._constraintZ4_cutoff = value

    @property
    def Splice02max_cutoff(self):
        return self._Splice02max_cutoff

    @Splice02max_cutoff.setter
    def Splice02max_cutoff(self, value):
        self._Splice02max_cutoff = value

    @property
    def Splice05max_cutoff(self):
        return self._Splice05max_cutoff

    @Splice05max_cutoff.setter
    def Splice05max_cutoff(self, value):
        self._Splice05max_cutoff = value

    @property
    def Splice08max_cutoff(self):
        return self._Splice08max_cutoff

    @Splice08max_cutoff.setter
    def Splice08max_cutoff(self, value):
        self._Splice08max_cutoff = value

    def categorize_variant(self, annotated_vcf: pd.DataFrame):
        result = defaultdict(int)

        for annotation_term_lists in self.annotate_each_variant(annotated_vcf):
            for combination in product(*annotation_term_lists):
                result[Category(*combination)] += 1

        return result

    def annotate_each_variant(self, annotated_vcf):
        """Newly annotated each variant using CWAS annotation terms.
        In order to annotate each variant with multiple annotation terms
        from each group efficiently, "Annotation integer" has defined.

        Annotation integer: A bitwise representation of the annotation
        of each variant where each bit means each annotation term

        e.g. If the annotation terms is ['A', 'B', 'C', 'D'] and
        the annotation integer is 0b1011, it means that
        the variant is annotated as 'A', 'B' and 'D'.
        """
        variant_type_annotation_ints = self.annotate_variant_type(annotated_vcf)
        conservation_annotation_ints = self.annotate_conservation(annotated_vcf)
        gene_list_annotation_ints = self.annotate_gene_list(annotated_vcf)
        gencode_annotation_ints = self.annotate_gencode(annotated_vcf)
        region_annotation_ints = self.annotate_region(annotated_vcf)

        for (
            variant_type_annotation_int,
            conservation_annotation_int,
            gene_list_annotation_int,
            gencode_annotation_int,
            region_annotation_int,
        ) in zip(
            variant_type_annotation_ints,
            conservation_annotation_ints,
            gene_list_annotation_ints,
            gencode_annotation_ints,
            region_annotation_ints,
        ):
            yield (
                self.parse_annotation_int(
                    variant_type_annotation_int, "variant_type",
                ),
                self.parse_annotation_int(
                    gene_list_annotation_int, "gene_list"
                ),
                self.parse_annotation_int(
                    conservation_annotation_int, "conservation",
                ),
                self.parse_annotation_int(gencode_annotation_int, "gencode"),
                self.parse_annotation_int(region_annotation_int, "region"),
            )

    def annotate_variant_type(self, annotated_vcf: pd.DataFrame) -> list:
        variant_type_annotation_idx = get_idx_dict(
            self._category_domain["variant_type"]
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

    def annotate_conservation(self, annotated_vcf: pd.DataFrame) -> list:
        conservation_annotation_idx = get_idx_dict(
            self._category_domain["conservation"]
        )
        phylop_conv_func = (
            lambda x: self._default_phylop
            if x == ""
            else max(map(float, x.split("&")))
        )
        phast_conv_func = (
            lambda x: self._default_phastcons
            if x == ""
            else max(map(float, x.split("&")))
        )
        jarvis099_conv_func = (
            lambda x: self._default_jarvis099
            if x == ""
            else max(map(float, x.split("&")))
        )
        ConstraintZ4_conv_func = (
            lambda x: self._default_constraintZ4
            if x == ""
            else max(map(float, x.split("&")))
        )
        Splice02max_conv_func = (
            lambda x: self._default_Splice02max
            if x == ""
            else max(map(float, x.split("&")))
        )
        Splice05max_conv_func = (
            lambda x: self._default_Splice05max
            if x == ""
            else max(map(float, x.split("&")))
        )
        Splice08max_conv_func = (
            lambda x: self._default_Splice08max
            if x == ""
            else max(map(float, x.split("&")))
        )

        phylop_scores = np.vectorize(phylop_conv_func)(
            annotated_vcf["phyloP46way"].values
        )
        phast_scores = np.vectorize(phast_conv_func)(
            annotated_vcf["phastCons46way"].values
        )
        jarvis099_scores = np.vectorize(jarvis099_conv_func)(
            annotated_vcf["JARVIS099"].values
        )
        ConstraintZ4_scores = np.vectorize(ConstraintZ4_conv_func)(
            annotated_vcf["ConstraintZ4"].values
        )
        Splice02max_scores = np.vectorize(Splice02max_conv_func)(
            annotated_vcf["Splice02max"].values
        )
        Splice05max_scores = np.vectorize(Splice05max_conv_func)(
            annotated_vcf["Splice05max"].values
        )
        Splice08max_scores = np.vectorize(Splice08max_conv_func)(
            annotated_vcf["Splice08max"].values
        )

        is_phylop_conservation_arr = (
            phylop_scores >= self._phylop_cutoff
        ).astype(np.int32)
        is_phast_conservation_arr = (
            phast_scores >= self._phastcons_cutoff
        ).astype(np.int32)
        is_jarvis099_conservation_arr = (
            jarvis099_scores >= self._jarvis099_cutoff
        ).astype(np.int32)
        is_ConstraintZ4_conservation_arr = (
            ConstraintZ4_scores >= self._constraintZ4_cutoff
        ).astype(np.int32)
        is_Splice02max_conservation_arr = (
            Splice02max_scores >= self._Splice02max_cutoff
        ).astype(np.int32)
        is_Splice05max_conservation_arr = (
            Splice05max_scores >= self._Splice05max_cutoff
        ).astype(np.int32)
        is_Splice08max_conservation_arr = (
            Splice08max_scores >= self._Splice08max_cutoff
        ).astype(np.int32)

        annotation_ints = np.vectorize(
            lambda is_phylop_conservation: 2
            ** conservation_annotation_idx["phyloP46way"]
            if is_phylop_conservation
            else 0
        )(is_phylop_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_phast_conservation: 2
            ** conservation_annotation_idx["phastCons46way"]
            if is_phast_conservation
            else 0
        )(is_phast_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_jarvis099_conservation: 2
            ** conservation_annotation_idx["JARVIS099"]
            if is_jarvis099_conservation
            else 0
        )(is_jarvis099_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_ConstraintZ4_conservation: 2
            ** conservation_annotation_idx["ConstraintZ4"]
            if is_ConstraintZ4_conservation
            else 0
        )(is_ConstraintZ4_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_Splice02max_conservation: 2
            ** conservation_annotation_idx["Splice02max"]
            if is_Splice02max_conservation
            else 0
        )(is_Splice02max_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_Splice05max_conservation: 2
            ** conservation_annotation_idx["Splice05max"]
            if is_Splice05max_conservation
            else 0
        )(is_Splice05max_conservation_arr)
        annotation_ints += np.vectorize(
            lambda is_Splice08max_conservation: 2
            ** conservation_annotation_idx["Splice08max"]
            if is_Splice08max_conservation
            else 0
        )(is_Splice08max_conservation_arr)
        annotation_ints += 2 ** conservation_annotation_idx["All"]

        return annotation_ints

    def annotate_gene_list(self, annotated_vcf: pd.DataFrame) -> list:
        gene_list_annotation_idx = get_idx_dict(
            self._category_domain["gene_list"]
        )
        gene_symbols = annotated_vcf["SYMBOL"].values
        gene_nearests = annotated_vcf["NEAREST"].values
        gencodes = annotated_vcf["Consequence"].values  # GENCODE annotations

        annotation_int_list = []
        annotation_int_dict = {}

        for symbol, nearest, gencode in zip(
            gene_symbols, gene_nearests, gencodes
        ):
            gene = (
                nearest
                if "downstream_gene_variant" in gencode
                or "intergenic_variant" in gencode
                else symbol
            )
            annotation_int = annotation_int_dict.get(gene, 0)

            if annotation_int == 0:
                gene_list_set = self._gene_matrix.get(gene, set())

                if gene_list_set:
                    for gene_cat in gene_list_annotation_idx:
                        if gene_cat in gene_list_set:
                            annotation_int += (
                                2 ** gene_list_annotation_idx[gene_cat]
                            )

            annotation_int_list.append(annotation_int)

        annotation_ints = np.asarray(annotation_int_list)
        annotation_ints += 2 ** gene_list_annotation_idx["Any"]

        return annotation_ints

    def annotate_gencode(self, annotated_vcf: pd.DataFrame) -> list:
        gencode_annotation_idx = get_idx_dict(self._category_domain["gencode"])
        gene_symbols = annotated_vcf["SYMBOL"].values
        gene_nearests = annotated_vcf["NEAREST"].values
        gencodes = annotated_vcf["Consequence"].values
        polyphens = annotated_vcf["PolyPhen"].values

        annotation_int_list = []

        for symbol, nearest, gencode, polyphen in zip(
            gene_symbols, gene_nearests, gencodes, polyphens
        ):
            gene = (
                nearest
                if "downstream_gene_variant" in gencode
                or "intergenic_variant" in gencode
                else symbol
            )
            gene_list_set = self._gene_matrix.get(gene, set())
            annotation_int = 0
            is_in_coding = False

            if "ProteinCoding" in gene_list_set:
                is_in_coding = True
                annotation_int += 2 ** gencode_annotation_idx["CodingRegion"]

                # Coding region
                if (
                    "stop_gained" in gencode
                    or "splice_donor" in gencode
                    or "splice_acceptor" in gencode
                ):
                    annotation_int += 2 ** gencode_annotation_idx["LoFRegion"]
                elif (
                    "frameshift_variant" in gencode
                    or "transcript_amplification" in gencode
                    or "transcript_ablation" in gencode
                ):
                    annotation_int += 2 ** gencode_annotation_idx["LoFRegion"]
                    annotation_int += (
                        2 ** gencode_annotation_idx["FrameshiftRegion"]
                    )
                elif (
                    "missense_variant" in gencode
                    or "start_lost" in gencode
                    or "stop_lost" in gencode
                ):
                    annotation_int += (
                        2 ** gencode_annotation_idx["MissenseRegion"]
                    )

                    if "probably_damaging" in polyphen:
                        annotation_int += (
                            2
                            ** gencode_annotation_idx[
                                "MissenseHVARDRegionSimple"
                            ]
                        )

                elif (
                    "inframe_deletion" in gencode
                    or "inframe_insertion" in gencode
                ):
                    annotation_int += (
                        2 ** gencode_annotation_idx["InFrameRegion"]
                    )
                elif "synonymous_variant" in gencode:
                    annotation_int += (
                        2 ** gencode_annotation_idx["SilentRegion"]
                    )
                elif (
                    "stop_retained_variant" not in gencode
                    and "incomplete_terminal_codon_variant" not in gencode
                    and "protein_altering_variant" not in gencode
                    and "coding_sequence_variant" not in gencode
                ):
                    # Noncoding
                    annotation_int = 0
                    is_in_coding = False

            if not is_in_coding:
                annotation_int += 2 ** gencode_annotation_idx["NoncodingRegion"]

                if "5_prime_UTR_variant" in gencode:
                    annotation_int += 2 ** gencode_annotation_idx["5PrimeUTRsRegion"]
                elif "3_prime_UTR_variant" in gencode:
                    annotation_int += 2 ** gencode_annotation_idx["3PrimeUTRsRegion"]
                elif "upstream_gene_variant" in gencode:
                    annotation_int += (
                        2 ** gencode_annotation_idx["PromoterRegion"]
                    )
                elif "intron_variant" in gencode:
                    annotation_int += (
                        2 ** gencode_annotation_idx["IntronRegion"]
                    )
                elif "splice_region_variant" in gencode:
                    annotation_int += (
                        2 ** gencode_annotation_idx["SpliceSiteNoncanonRegion"]
                    )
                elif (
                    "downstream_gene_variant" in gencode
                    or "intergenic_variant" in gencode
                ):
                    annotation_int += (
                        2 ** gencode_annotation_idx["IntergenicRegion"]
                    )
                elif "ProteinCoding" not in gene_list_set:
                    if "Antisense" in gene_list_set:
                        annotation_int += (
                            2 ** gencode_annotation_idx["AntisenseRegion"]
                        )
                    elif "lincRNA" in gene_list_set:
                        annotation_int += (
                            2 ** gencode_annotation_idx["lincRnaRegion"]
                        )
                    else:
                        annotation_int += (
                            2 ** gencode_annotation_idx["OtherTranscriptRegion"]
                        )

            annotation_int_list.append(annotation_int)

        annotation_ints = np.asarray(annotation_int_list)
        annotation_ints += 2 ** gencode_annotation_idx["Any"]

        return annotation_ints

    def annotate_region(self, annotated_vcf: pd.DataFrame) -> list:
        region_annotation_idx = get_idx_dict(self._category_domain["region"])
        annotation_ints = np.zeros(len(annotated_vcf.index), dtype=int)

        for region in region_annotation_idx:
            if region == "Any":
                continue

            region_vals = annotated_vcf[region].values.astype(np.int32)
            annotation_int_conv_func = (
                lambda x: 2 ** region_annotation_idx[region] * x
            )
            annotation_ints += np.vectorize(annotation_int_conv_func)(
                region_vals
            )

        annotation_ints += 2 ** region_annotation_idx["Any"]

        return annotation_ints

    def parse_annotation_int(
        self, annotation_int: int, annotation_term_type: str
    ) -> list:
        """ Parse the annotation integer and
        choose the appropriate subset from the specific annotation terms.
        """
        return extract_sublist_by_int(
            self._category_domain[annotation_term_type], annotation_int
        )

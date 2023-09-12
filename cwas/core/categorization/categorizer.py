"""
Module for the categorization of the variants into CWAS categories and
counting the number of variants in each category.
CWAS category is a combination of annotation terms
from each of annotation groups.

There are currently 5 groups of annotation terms.

--- The groups of the annotation terms for CWAS ---
    1. Variant types (variant_type)
    2. Functional score (functional_score)
    3. Gene lists (gene_set)
    4. GENCODE annotation categories (gencode)
    5. Functional annotation categories (functional_annotation)

"""
from collections import defaultdict
from itertools import product

import numpy as np
import pandas as pd
from tqdm import tqdm
from cwas.core.categorization.utils import extract_sublist_by_int, get_idx_dict


class Categorizer:
    def __init__(self, category_domain: dict, gene_matrix: dict, mis_info_key: str, mis_thres: float) -> None:
        self._category_domain = category_domain
        self._gene_matrix = gene_matrix
        self._mis_info_key = mis_info_key
        self._mis_thres = mis_thres

    def categorize_variant(self, annotated_vcf: pd.DataFrame):
        result = defaultdict(int)

        for annotation_term_lists in self.annotate_each_variant(annotated_vcf):
            for combination in product(*annotation_term_lists):
                result["_".join(combination)] += 1

        return result

    def get_intersection(self, annotated_vcf: pd.DataFrame):
        result = defaultdict(lambda: defaultdict(int))

        for annotation_term_lists in tqdm(self.annotate_each_variant(annotated_vcf), total=len(annotated_vcf)):
            for combination in product(*annotation_term_lists):
                for combination2 in product(*annotation_term_lists):
                    result["_".join(combination)]["_".join(combination2)] += 1

        return result

    def get_intersection_variant_level(self, annotated_vcf, category_combinations):
      # Generate unique category combinations
      #category_combinations = set()
      #for annotation_term_lists in annotate_each_variant(annotated_vcf):
      #    for combination in product(*annotation_term_lists):
      #        category_combinations.add("_".join(combination))
      
      # Create a matrix with zeros
      num_variants = annotated_vcf.shape[0]
      num_categories = len(category_combinations)
      #matrix = np.zeros((num_variants, num_categories), dtype=int)
      matrix = np.zeros((num_variants, num_categories), dtype='uint64')
      
      # Create a dictionary to map categories to matrix column indices
      category_to_index = {category: index for index, category in enumerate(category_combinations)}

      # Populate the matrix
      for variant_index, annotation_term_lists in enumerate(self.annotate_each_variant(annotated_vcf)):
          for combination in product(*annotation_term_lists):
              category = "_".join(combination)
              if category in category_to_index:
                  category_index = category_to_index[category]
                  matrix[variant_index, category_index] = 1

      df = pd.DataFrame(matrix, columns=category_combinations)

      return df

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
        #functional_score_annotation_ints = self.annotate_functional_score(annotated_vcf)
        gene_set_annotation_ints = self.annotate_gene_set(annotated_vcf)
        gencode_annotation_ints = self.annotate_gencode(annotated_vcf)
        #region_annotation_ints = self.annotate_region(annotated_vcf)
        functional_ints = annotated_vcf['ANNOT']

        for (
            variant_type_annotation_int,
            functional_score_annotation_int,
            gene_set_annotation_int,
            gencode_annotation_int,
            region_annotation_int,
        ) in zip(
            variant_type_annotation_ints,
            functional_ints,
            gene_set_annotation_ints,
            gencode_annotation_ints,
            functional_ints,
        ):
            yield (
                self.parse_annotation_int(
                    variant_type_annotation_int, "variant_type",
                ),
                self.parse_annotation_int(
                    gene_set_annotation_int, "gene_set"
                ),
                self.parse_annotation_int_(
                    functional_score_annotation_int, "functional_score",
                ),
                self.parse_annotation_int(gencode_annotation_int, "gencode"),
                self.parse_annotation_int_(region_annotation_int, "functional_annotation"),
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

    def annotate_functional_score(self, annotated_vcf: pd.DataFrame) -> list:
        functional_score_annotation_idx = get_idx_dict(
            self._category_domain["functional_score"]
        )
        annotation_ints = np.zeros(len(annotated_vcf.index), dtype=int)

        for score in functional_score_annotation_idx:
            if score == "All":
                continue

            score_vals = annotated_vcf[score].values.astype(np.int32)
            annotation_int_conv_func = (
                lambda x: 2 ** functional_score_annotation_idx[score] * x
            )
            annotation_ints += np.vectorize(annotation_int_conv_func)(
                score_vals
            )

        annotation_ints += 2 ** functional_score_annotation_idx["All"]

        return annotation_ints

    def annotate_gene_set(self, annotated_vcf: pd.DataFrame) -> list:
        gene_set_annotation_idx = get_idx_dict(
            self._category_domain["gene_set"]
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
                gene_set_ = self._gene_matrix.get(gene, set())

                if gene_set_:
                    for gene_cat in gene_set_annotation_idx:
                        if gene_cat in gene_set_:
                            annotation_int += (
                                2 ** gene_set_annotation_idx[gene_cat]
                            )

            annotation_int_list.append(annotation_int)

        annotation_ints = np.asarray(annotation_int_list)
        annotation_ints += 2 ** gene_set_annotation_idx["Any"]

        return annotation_ints

    def annotate_gencode(self, annotated_vcf: pd.DataFrame) -> list:
        gencode_annotation_idx = get_idx_dict(self._category_domain["gencode"])
        gene_symbols = annotated_vcf["SYMBOL"].values
        gene_nearests = annotated_vcf["NEAREST"].values
        gencodes = annotated_vcf["Consequence"].values
        lofs = annotated_vcf["LoF"].values
        lof_flags = annotated_vcf["LoF_flags"].values
        mis_scores = annotated_vcf["MisDb_" + self._mis_info_key].values

        annotation_int_list = []

        for symbol, nearest, gencode, lof, lof_flag, mis_score in zip(
            gene_symbols, gene_nearests, gencodes, lofs, lof_flags, mis_scores
        ):
            gene = (
                nearest
                if "downstream_gene_variant" in gencode
                or "intergenic_variant" in gencode
                else symbol
            )
            gene_set_ = self._gene_matrix.get(gene, set())
            annotation_int = 0
            is_in_coding = False

            if "ProteinCoding" in gene_set_:
                is_in_coding = True
                annotation_int += 2 ** gencode_annotation_idx["CodingRegion"]

                # Coding region
                if ((
                    "stop_gained" in gencode
                    or "splice_donor" in gencode
                    or "splice_acceptor" in gencode
                )
                and (lof == "HC")
                and ((lof_flag=='SINGLE_EXON') or (lof_flag==""))
                ):
                    annotation_int += 2 ** gencode_annotation_idx["PTVRegion"]
                elif ((
                    "frameshift_variant" in gencode
                )
                and (lof == "HC")
                and ((lof_flag=='SINGLE_EXON') or (lof_flag==""))
                ):
                    annotation_int += 2 ** gencode_annotation_idx["PTVRegion"]
                    annotation_int += (
                        2 ** gencode_annotation_idx["FrameshiftRegion"]
                    )
                elif (
                    "missense_variant" in gencode
                    or "protein_altering_variant" in gencode
                    or "start_lost" in gencode
                    or "stop_lost" in gencode
                ):
                    annotation_int += (
                        2 ** gencode_annotation_idx["MissenseRegion"]
                    )

                    if ((len(mis_score)!=0)
                        and (mis_score!='NA')):
                        if '&' in mis_score:
                            mis_score = max([float(i) for i in mis_score.split('&')])
                        else:
                            mis_score = float(mis_score)
                        if mis_score >= self._mis_thres:
                            annotation_int += (
                                2 ** gencode_annotation_idx["DamagingMissenseRegion"]
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
                    and "stop_gained" not in gencode
                    and "splice_donor" not in gencode
                    and "splice_acceptor" not in gencode
                    and "frameshift_variant" not in gencode
                ):
                    # Noncoding
                    annotation_int = 0
                    is_in_coding = False

            if not is_in_coding:
                annotation_int += 2 ** gencode_annotation_idx["NoncodingRegion"]

                if "_UTR_" in gencode:
                    annotation_int += 2 ** gencode_annotation_idx["UTRsRegion"]
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
                elif "ProteinCoding" not in gene_set_:
                    if "lincRNA" in gene_set_:
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
        region_annotation_idx = get_idx_dict(self._category_domain["functional_annotation"])
        #annotation_floats = np.zeros(len(annotated_vcf.index), dtype=float)
        annotation_ints = np.zeros(len(annotated_vcf.index), dtype=int)

        for region in region_annotation_idx:
            if region == "Any":
                continue

            region_vals = annotated_vcf[region].values.astype(np.int32)
            annotation_int_conv_func = (
                lambda x: 2 ** int(region_annotation_idx[region]) * int(x)
            )
            annotation_ints = annotation_ints + np.vectorize(annotation_int_conv_func)(region_vals)
            #annotation_floats += np.vectorize(annotation_int_conv_func, otypes=[float])(
            #    region_vals
            #)

        #annotation_ints = np.array([int(x) for x in annotation_floats])
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

    def parse_annotation_int_(
        self, annotation_int: int, annotation_term_type: str
    ) -> list:
        """ Parse the annotation integer and
        choose the appropriate subset from the specific annotation terms.
        """
        cat_list = [*self._category_domain['functional_score'], *self._category_domain['functional_annotation']]
        cat_list = [item for item in cat_list if item not in ['Any', 'All']]
        labels = extract_sublist_by_int(
            cat_list, int(annotation_int)
        )
        if annotation_term_type == 'functional_score':
            return ['All'] + [item for item in labels if item in self._category_domain['functional_score']]
        elif annotation_term_type == 'functional_annotation':
            return ['Any'] + [item for item in labels if item in self._category_domain['functional_annotation']]
        
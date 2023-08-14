class Category:
    def __init__(
        self, variant_type, gene_set, functional_score, gencode, functional_annotation
    ) -> None:
        self._variant_type = variant_type
        self._gene_set = gene_set
        self._functional_score = functional_score
        self._gencode = gencode
        self._functional_annotation = functional_annotation

    def __repr__(self):
        return (
            f"Category("
            f'"{self._variant_type}", '
            f'"{self._gene_set}", '
            f'"{self._functional_score}", '
            f'"{self._gencode}", '
            f'"{self._functional_annotation}")'
        )

    def __str__(self):
        """This representation is from An et al., 2018"""
        return "_".join(
            [
                self._variant_type,
                self._gene_set,
                self._functional_score,
                self._gencode,
                self._functional_annotation,
            ]
        )

    def __eq__(self, other):
        return (
            self._variant_type == other._variant_type
            and self._gene_set == other._gene_set
            and self._functional_score == other._functional_score
            and self._gencode == other._gencode
            and self._functional_annotation == other._functional_annotation
        )

    def __hash__(self):
        return hash(repr(self))

    def to_dict(self):
        return {
            "variant_type": self._variant_type,
            "gene_set": self._gene_set,
            "functional_score": self._functional_score,
            "gencode": self._gencode,
            "functional_annotation": self._functional_annotation,
        }

    @staticmethod
    def from_str(category_str):
        return Category(*category_str.split("_"))

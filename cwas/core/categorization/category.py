class Category:
    def __init__(
        self, variant_type, gene_list, conservation, gencode, region
    ) -> None:
        self._variant_type = variant_type
        self._gene_list = gene_list
        self._conservation = conservation
        self._gencode = gencode
        self._region = region

    def __repr__(self):
        return (
            f"Category("
            f'"{self._variant_type}", '
            f'"{self._gene_list}", '
            f'"{self._conservation}", '
            f'"{self._gencode}", '
            f'"{self._region}")'
        )

    def __str__(self):
        """This representation is from An et al., 2018"""
        return "_".join(
            [
                self._variant_type,
                self._gene_list,
                self._conservation,
                self._gencode,
                self._region,
            ]
        )

    def __eq__(self, other):
        return (
            self._variant_type == other._variant_type
            and self._gene_list == other._gene_list
            and self._conservation == other._conservation
            and self._gencode == other._gencode
            and self._region == other._region
        )

    def __hash__(self):
        return hash(repr(self))

    def to_dict(self):
        return {
            "variant_type": self._variant_type,
            "gene_list": self._gene_list,
            "conservation": self._conservation,
            "gencode": self._gencode,
            "region": self._region,
        }

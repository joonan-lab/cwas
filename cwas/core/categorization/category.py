class Category:
    def __init__(
        self, variant_type, conservation, gene_list, gencode, region
    ) -> None:
        self._variant_type = variant_type
        self._conservation = conservation
        self._gene_list = gene_list
        self._gencode = gencode
        self._region = region

    def __repr__(self):
        return "&".join(
            [
                self._variant_type,
                self._conservation,
                self._gene_list,
                self._gencode,
                self._region,
            ]
        )

    def __eq__(self, other):
        return (
            self._variant_type == other._variant_type
            and self._conservation == other._conservation
            and self._gene_list == other._gene_list
            and self._gencode == other._gencode
            and self._region == other._region
        )

    def __hash__(self):
        return hash(repr(self))

    def to_dict(self):
        return {
            'variant_type': self._variant_type,
            'conservation': self._conservation,
            'gene_list': self._gene_list,
            'gencode': self._gencode,
            'region': self._region
        }

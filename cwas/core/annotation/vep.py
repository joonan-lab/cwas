from cwas.utils.check import check_is_file


class VEP:
    def __init__(self, vep: str) -> None:
        self.vep = vep
        self._check_vep_path()

    def _check_vep_path(self):
        check_is_file(self.vep)

    def get_vep_path(self):
        return self.vep

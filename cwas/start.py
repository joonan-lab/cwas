import argparse
import shutil
from pathlib import Path
from typing import Optional

import cwas.utils.log as log
from cwas.runnable import Runnable


class Start(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self.config_path = self.workspace / "configuration.txt"

    @staticmethod
    def _create_arg_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Arguments for Initializing a CWAS workspace",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        default_workspace = Path.home() / ".cwas"
        parser.add_argument(
            "-w",
            "--workspace",
            dest="workspace",
            required=False,
            type=Path,
            default=default_workspace,
            help="Path to your CWAS workspace directory",
        )
        return parser

    @staticmethod
    def _print_args(args: argparse.Namespace):
        log.print_arg("CWAS Workspace", args.workspace)

    @property
    def workspace(self):
        return self.args.workspace.resolve()

    def run(self):
        if self.workspace.is_dir():
            log.print_warn(
                f"The workspace '{self.workspace}' already exists. "
                f"This will be not re-created."
            )
        else:
            self._create_workspace()

        self._update_env()

        if self.config_path.is_file():
            log.print_warn(
                f"The configuration file '{self.config_path}' already exists. "
                f"This will be not re-created"
            )
        else:
            self._create_config_file()

    def _create_workspace(self):
        log.print_progress(f"Create CWAS workspace '{self.workspace}'")
        try:
            self.workspace.mkdir(exist_ok=True)
        except NotADirectoryError:
            log.print_err("The path to CWAS workspace is invalid.")
            raise

    def _update_env(self):
        self.env.set_env("CWAS_WORKSPACE", self.workspace)
        self.env.save()

    def _create_config_file(self):
        config = self._init_config()
        with self.config_path.open("w") as config_file:
            for k, v in config.items():
                print(f"{k}={v}", file=config_file)

    def _init_config(self) -> dict:
        config = {
            "ANNOTATION_DATA_DIR": "",
            "GENE_MATRIX": "",
            "ANNOTATION_KEY_CONFIG": "",
            "VEP": "",
            "VEP_CACHE_DIR": "",
            "VEP_CONSERVATION_FILE": "",
            "VEP_LOFTEE": "",
            "VEP_HUMAN_ANCESTOR_FA": "",
            "VEP_GERP_BIGWIG": "",
            "VEP_MIS_DB": "",
            "VEP_MIS_INFO_KEY": "",
            "VEP_MIS_THRES": "",
        }

        installed_vep = self._find_vep_path()
        if installed_vep:
            config["VEP"] = installed_vep

        return config

    def _find_vep_path(self) -> Optional[str]:
        log.print_progress("Find pre-installed VEP")
        vep = shutil.which("vep")
        log.print_progress(f"VEP path: '{vep}'")
        return vep
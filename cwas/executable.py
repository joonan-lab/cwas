import argparse


class Executable:
    def run(self, args: argparse.ArgumentParser):
        pass

    def create_arg_parser(self) -> argparse.ArgumentParser:
        pass

    def print_args(self, args: argparse.ArgumentParser):
        pass

    def check_args_validity(self, args: argparse.ArgumentParser):
        pass

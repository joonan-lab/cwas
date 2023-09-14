"""
Command line tool for Category-wide association study (CWAS)
"""
import sys

import cwas.factory
from cwas.utils.log import print_log


def main():
    # TODO: Add conditional statements or exception handling. This is very
    #  error-prone.
    print_log("LOG", "Category-Wide Association Study (CWAS)")
    cwas_factory = cwas.factory.create(sys.argv[1])
    cwas_args = cwas_factory.argparser().parse_args(sys.argv[2:])
    cwas_obj = cwas_factory.runnable
    print_log("LOG", f"Current step: {cwas_obj.__name__}")
    cwas_inst = cwas_obj(cwas_args)
    return cwas_inst
    #cwas_inst.run()

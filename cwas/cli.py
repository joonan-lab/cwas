"""
Command line tool for CWAS
"""
import sys

from cwas.factory import get_runnable


def main():
    # TODO: Add conditional statements or exception handling. This is very
    #  error-prone.
    cwas_step_name = sys.argv[1]
    cwas_args = sys.argv[2:]
    cwas_obj = get_runnable(cwas_step_name)
    cwas_inst = cwas_obj.get_instance(cwas_args)
    cwas_inst.run()

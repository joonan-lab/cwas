"""
Command line tool for Category-wide association study (CWAS)
"""
import sys

from cwas.env import Env
from cwas.factory import get_runnable
from cwas.utils.log import print_log


def main():
    # TODO: Add conditional statements or exception handling. This is very
    #  error-prone.
    print_log("LOG", "Category-Wide Association Study (CWAS)")
    cwas_env = Env()
    cwas_step_name = sys.argv[1]
    cwas_args = sys.argv[2:]
    cwas_obj = get_runnable(cwas_step_name)
    print_log("LOG", f"Current step: {cwas_obj.__name__}")
    cwas_inst = cwas_obj.get_instance(cwas_args)
    setattr(cwas_inst, "env", cwas_env)
    cwas_inst.run()
    cwas_env.save()

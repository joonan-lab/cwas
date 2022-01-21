import random
from multiprocessing import cpu_count

import pytest
from cwas.categorization import Categorization


def test_parse_args():
    cpu = random.choice(range(1, cpu_count() + 1))
    inst = Categorization.get_instance(["-p", str(cpu)])
    assert inst.num_proc == cpu


def test_parse_args_invalid_num_proc():
    with pytest.raises(ValueError):
        Categorization.get_instance(["-p", str(cpu_count() + 1)])

    with pytest.raises(ValueError):
        Categorization.get_instance(["-p", "0"])

    with pytest.raises(ValueError):
        Categorization.get_instance(["-p", "-1"])


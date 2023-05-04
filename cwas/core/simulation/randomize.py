import numpy as np

def label_variant(ref: str, alt: str) -> int:
    """ Return an integer according to the type of the input small variant

    ** Labels **
        INDEL1 if a difference of lengths between ref and alt is 3n + 1
        INDEL2 if the difference is 3n + 2
        INDEL3 if the difference is 3n + 3
        where n is 0 or a positive integer
    """
    assert len(ref) == 1 or len(alt) == 1, \
        'This function current support only small variants ' \
        'such as SNV and INDEL.'
    len_diff = abs(len(ref) - len(alt))

    if len_diff == 0:
        return 0  # SNV
    elif len_diff % 3 == 1:
        return 1  # INDEL1
    elif len_diff % 3 == 2:
        return 2  # INDEL2
    else:  # len_diff % 3 == 0
        return 3  # INDEL3
    
def pick_mutation() -> (str, str):
    """ Get a mutation from the mutation distribution model (Lynch 2010). """
    
    x = np.random.uniform(0, 1)

    if x <= 0.211062887:
        ref = 'C'
        alt = 'T'
    elif x <= 0.422125774:
        ref = 'G'
        alt = 'A'
    elif x <= 0.551200326:
        ref = 'A'
        alt = 'G'
    elif x <= 0.680274877:
        ref = 'T'
        alt = 'C'
    elif x <= 0.728393387:
        ref = 'G'
        alt = 'T'
    elif x <= 0.776511898:
        ref = 'C'
        alt = 'A'
    elif x <= 0.821985623:
        ref = 'G'
        alt = 'C'
    elif x <= 0.867459349:
        ref = 'C'
        alt = 'G'
    elif x <= 0.900590744:
        ref = 'T'
        alt = 'A'
    elif x <= 0.933722139:
        ref = 'A'
        alt = 'T'
    elif x <= 0.96686107:
        ref = 'A'
        alt = 'C'
    else:
        ref = 'T'
        alt = 'G'

    return ref, alt
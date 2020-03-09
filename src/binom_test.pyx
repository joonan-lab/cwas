import numpy as np
import scipy.stats as ss


cpdef binomTest(np.ndarray[long] vals1, np.ndarray[long] vals2):
    cdef int n1, n2
    out = []
    for n1, n2 in zip(vals1, vals2):
        out.append(ss.binom_test(x=n1, n=n1 + n2, p=0.5, alternative='two-sided'))
    return out


cpdef binomTest_onesided(np.ndarray[long] vals1, np.ndarray[long] vals2):
    cdef int n1, n2
    out = []
    for n1, n2 in zip(vals1, vals2):
        if n1 > n2:
            out.append(ss.binom_test(x=n1, n=n1 + n2, p=0.5, alternative='greater'))
        else:
            out.append(ss.binom_test(x=n2, n=n1 + n2, p=0.5, alternative='greater'))
    return out


cpdef check_pToZ(np.ndarray[double] p):
    cdef np.ndarray[double] p1
    cdef np.ndarray[double] z
    # Transform to 1 - p
    p1 = 1 - p
    print(p1[:10])
    # Convert to z score
    z = ss.norm.ppf(p1)
    print(z[:10])
    return (z)

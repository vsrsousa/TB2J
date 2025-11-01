import sys

import numpy as np

MAX_EXP_ARGUMENT = np.log(sys.float_info.max / 100000)
print(MAX_EXP_ARGUMENT)


def fermi(e, mu, width=0.01, nspin=1):
    """
    the fermi function.
     .. math::
        f=\\frac{nspin}{\exp((e-\mu)/width)+1}

    :param e,mu,width: e,\mu,width
    :param nspin: spin degeneracy factor (1 for spin-polarized, 2 for non-spin-polarized)
    """
    x = (e - mu) / width
    # disable overflow warning
    # with warnings.catch_warnings():
    #    warnings.simplefilter("ignore")
    #    ret = np.where(x < MAX_EXP_ARGUMENT, 1 / (1.0 + np.exp(x)), 0.0)

    ret = np.zeros_like(x, dtype=float)
    for i, xi in enumerate(x):
        if xi < 700:
            ret[i] = 1 / (1.0 + np.exp(xi))
    return ret * nspin

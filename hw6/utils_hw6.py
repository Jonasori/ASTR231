"""
Some helper functions for hw6.

blah
"""


# Get a chi2 value
def chi2(data, model):
    # Takes two lists.
    c = 0
    for i in range(len(data)):
        c += (data[i]-model[i])**2
    return c


def round(num, n_decs):
    i = 0
    while i < len(str(num)):
        if str(num)[i] == '.':
            rounded_num_str = str(num)[:i] + str(num)[:i + 1 + n_decs]
            return float(rounded_num_str)
        i += 1

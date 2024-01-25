import matplotlib.pyplot as plt
import numpy as np

"""
Test implications of using positive or negative root in solution of the 1D outflow problem (particularly BCs).
"""


def calc_n(z, delta, neg_root=False):
    c = delta + 1 / delta
    R = np.sqrt((delta + 1 / delta) ** 2 - 4 * z**2)
    if neg_root:
        R *= -1
    return (c + R) / 2


def main(delta=0.1, zmin=-1, zmax=1, nz=200):
    z = np.linspace(zmin, zmax, nz)
    npos = calc_n(z, delta)
    nneg = calc_n(z, delta, neg_root=True)

    plt.plot(z, npos, color="r", label="+ve root")
    plt.plot(z, nneg, color="b", linestyle="--", label="-ve root")
    plt.legend()
    plt.show()


main(delta=0.7)

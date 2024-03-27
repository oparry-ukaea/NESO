import numpy as np
import matplotlib.pyplot as plt
import random


def gen_grid(ranges=[[0, 16e-3], [0, 16e-3], [0, 10]], Ncells=[64, 64, 64], dims=3):
    lin_pts = {}
    lbls = ["x", "y", "z"][0:dims]
    for r, N, lbl in zip(ranges, Ncells, lbls):
        lin_pts[lbl] = np.linspace(r[0], r[1], N)
    if dims == 3:
        return np.meshgrid(lin_pts["x"], lin_pts["y"], lin_pts["z"])
    elif dims == 2:
        return np.meshgrid(lin_pts["x"], lin_pts["y"])
    else:
        raise ValueError("gen_grid: Only set up for dims=2 or dims=3")


def gen_phases(nmodes):
    random.seed(1)
    return [random.uniform(-np.pi, np.pi) for imode in range(nmodes)]


def mixmode(x, nmodes=14, peak_mode=4):
    phases = gen_phases(nmodes)
    result = np.zeros_like(x)
    for imode, phase in enumerate(phases):
        # N.B. BOUT++ docs say 1-14, but code does 1-13
        prefac = 1.0 / (1.0 + np.abs(imode - peak_mode))
        result += prefac * prefac * np.cos(imode * x + phase)
    return result


def test_w(x, y, z=1):
    twoPI = 2 * np.pi
    Nx = 10
    Ny = 3
    xs = x / 16e-3 * Nx * twoPI
    ys = y / 16e-3 * Ny * twoPI
    return np.sin(xs) * np.sin(ys)


def main():
    x, y = gen_grid(dims=2)
    z = np.zeros_like(x)
    w1 = mixmode(x * 2 * np.pi / 16e-3, peak_mode=4)
    w2 = mixmode(y / 16e-3 - z, peak_mode=2)
    w = w1 * w2
    plt.pcolor(x, y, w)
    plt.show(block=True)


main()

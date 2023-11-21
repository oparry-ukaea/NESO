import matplotlib.pyplot as plt
import numpy as np

mu = 1
sigma = 0.1
prefac = 0.8
offset = 1.0
xmin = 0.0
xmax = 2.0
x = np.arange(0.0, 2.0, 1e-2)
n_pert = offset + prefac * np.exp(-(x - mu) * (x - mu) / sigma / sigma)


delta = 0.9
n = 0.5 * ((delta + 1 / delta) + np.sqrt((delta + 1 / delta) ** 2 - 4 * (x - 1) ** 2))
G = x - 1
u = G / n
# print(f"{n_pert[0]:.18f}, {n_pert[1]:.18f}")
# plt.plot(x, n_pert)
# plt.plot(x, n, linestyle="--")
# plt.hlines([delta + 1 / delta], [xmin], [xmax])

print(f"max vel is {np.max(u):.3f}")
# print(delta + 1 / delta)

# delta_arr = np.arange(1e-2, 1.0, 1e-2)
# plt.plot(delta_arr, delta_arr + 1 / delta_arr)

plt.plot(x, u)

plt.show()

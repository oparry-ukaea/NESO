import math
import numpy as np


def get_norms(params, num_nucleons=2, verbose=False):

    # Constants
    Mp = 1.6726219e-27  # Proton mass [kg]
    Me = 9.1093837e-31  # Electron mass [kg]
    qe = 1.60217662e-19  # electron charge [C] or [J ev^-1]
    e0 = 8.854187817e-12  # Vacuum permittivity [F m^-1]
    kb = 1.38064852e-23  # Boltzmann constant [J K^-1]

    # Derived facs
    Mi = Mp * num_nucleons  # Ion mass  [kg]
    omega_ci = qe * params["B_T"] / Mi  # ion gyrofrequency [s^-1]
    clog = coulomb_log(params["T0_eV"], params["n0"], Me, Mi, params["Z"])
    rho_s0 = np.sqrt(Mi * params["T0_eV"] / (qe * params["B_T"]))

    if verbose:
        print("  Norm factors:")
        print(f"    rho_s0   = {rho_s0:.3e}")
        print(f"    omega_ci = {omega_ci:.3e}")

    eta = 5.2e-5 * params["Z"] * clog / params["T0_eV"] ** 1.5
    alpha = params["T0_eV"] / params["n0"] / qe / eta / omega_ci
    k = 2 * math.pi / params["Dz"]
    alpha2D = alpha * k**2
    print("\n  HW params:")
    print(f"    alpha   = {alpha:.2e}")
    print(f"    alpha2D = {alpha2D:.2e}")
    kappa = rho_s0 / params["lambda_q"]
    print(f"    kappa   = {kappa:.2e}")

    # Return norm facs
    norms = dict(dist=1 / rho_s0, time=omega_ci)
    return norms


def coulomb_log(T0, n0, Me, Mi, Z):
    if T0 < 0.1 or n0 < 1e10:
        clog = 10
    elif T0 * Me / Mi < T0 and T0 < 10 * Z**2:
        clog = 30 - 0.5 * np.log(n0) + 1.5 * np.log(T0)
    elif T0 * Me / Mi < 10 * Z**2 and 10 * Z**2 < T0:
        clog = 31 - 0.5 * np.log(n0) + np.log(T0)
    elif T0 < T0 * Me / Mi:
        clog = 23 - 0.5 * np.log(n0) + 1.5 * np.log(T0) - np.log(Z * Mi)
    return clog


def report_norm_params(params, fnorm):
    Dx = params["Dx"] * fnorm["dist"]
    Dy = params["Dy"] * fnorm["dist"]
    Dz = params["Dz"] * fnorm["dist"]
    dx = Dx / params["Nx"]
    dy = Dy / params["Ny"]
    dz = Dz / params["Nz"]
    print("\n  Normalised quantities: ")
    print(f"     Mesh dims: [{Dx:.2e}, {Dy:.2e}, {Dz:.2e}]")
    print(f"     Cell dims: [{dx:.2e}, {dy:.2e}, {dz:.2e}]")
    tTot_norm = params["tTot"] * fnorm["time"]
    print(f"    Total time: {tTot_norm:.2e}")
    dt_norm = tTot_norm / params["nsteps"]
    print(f"            dt: {dt_norm:.2e}")


def make_params(**kwargs):
    common_params = dict(
        Dx=18.2e-3,
        Dy=18.2e-3,
        Dz=10.0,
        Nx=64,
        Ny=64,
        Nz=64,
        T0_eV=40.0,
        B_T=0.5,
        lambda_q=6.462e-3,
        n0=1.59e19,
        Z=1,
        tTot=500 * 1e-6,
        nsteps=50000,
    )
    params = dict(common_params)
    params.update(**kwargs)
    return params


def main():
    # Params - dimensionless or SI

    all_params = dict(
        Case1=make_params(),
        Case2=make_params(n0=1.45e18),
        Case3=make_params(n0=1.33e17),
    )
    for lbl, params in all_params.items():
        print(f"\n{lbl}:")
        fnorm = get_norms(params, verbose=(lbl == "Case1"))
        report_norm_params(params, fnorm)


main()

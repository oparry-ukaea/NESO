
#include <boost/core/ignore_unused.hpp>

#include "LaxFriedrichsWithOutflowSolver.h"

namespace NESO::Solvers::H3LAPD {

std::string LaxFriedrichsWithOutflowSolver::solver_name =
    SU::GetRiemannSolverFactory().RegisterCreatorFunction(
        "LFOutflow", LaxFriedrichsWithOutflowSolver::create,
        "Lax-Friedrichs solver customised for outflows");

/**
 * @class LaxFriedrichsWithOutflowSolver
 *
 * @brief Customised Lax-Friedrichs Riemann solver.
 *
 */

/**
 * @brief Default constructor.
 */
LaxFriedrichsWithOutflowSolver::LaxFriedrichsWithOutflowSolver(
    const LU::SessionReaderSharedPtr &session)
    : SU::RiemannSolver(session) {}

/**
 *
 * @param Fwd   Forwards trace space.
 * @param Bwd   Backwards trace space.
 * @param flux  Resulting flux.
 */
void LaxFriedrichsWithOutflowSolver::v_Solve(const int nDim,
                                             const Nek2DArrConstInner &Fwd,
                                             const Nek2DArrConstInner &Bwd,
                                             Nek2DArr &flux) {

  ASSERTL1(CheckScalars("norms"), "Vn not defined.");
  const Nek1DArr &norms = m_scalars["norms"]();

  ASSERTL1(CheckScalars("z"), "z not defined.");
  const Nek1DArr &z = m_scalars["z"]();
  for (int j = 0; j < Fwd[0].size(); ++j) {
    // Field variables
    NekDouble n_L = Fwd[0][j];
    NekDouble nu_L = Fwd[1][j];
    NekDouble n_R = Bwd[0][j];
    NekDouble nu_R = Bwd[1][j];

    bool is_left = false;
    bool is_right = false;
    switch (nDim) {
    case 1:
      is_left = j == 0;
      is_right = j == 1;
      break;
    case 2:
      std::cout << "2D outflow not setup." << std::endl;
      break;
    case 3:
      const int nskip = 4 * 10 * 200 * 16;
      const int ntrace_outflow = 10 * 10 * 16;

      if (j >= nskip && j < nskip + ntrace_outflow) {
        is_left = true;
        if (z[j] > 1e-6 || norms[j] > -0.9999) {
          std::cout << j << ": z=" << z[j] << " with norm " << norms[j]
                    << " was assigned to LEFT boundary!?" << std::endl;
        }
      } else if (j >= nskip + ntrace_outflow &&
                 j < nskip + 2 * ntrace_outflow) {
        is_right = true;
        if (z[j] < 1.9999 || norms[j] < 0.9999) {
          std::cout << j << ": z=" << z[j] << " with norm " << norms[j]
                    << " was assigned to LEFT boundary!?" << std::endl;
        }
      }
      break;
    }

    if (is_left) {
      // x = 0
      n_R = n_L;
      // nu_R = nu_L;

      NekDouble u_L = nu_L / n_L;
      const NekDouble uinf = -1.0;
      NekDouble ub = u_L + (uinf - u_L * uinf * (-1.0)) * (-1.0);
      nu_R = n_L * (ub);
    } else if (is_right) {
      // x = 2
      n_R = n_L;
      // nu_R = nu_L;
      NekDouble u_L = nu_L / n_L;
      const NekDouble uinf = 1.0;
      NekDouble ub = u_L + (uinf - u_L * uinf * (1.0)) * (1.0);
      nu_R = n_L * (ub);
    }

    // Velocity
    NekDouble u_L = nu_L / n_L;
    NekDouble u_R = nu_R / n_R;

    // const temperature
    NekDouble T = 1.0;

    // sound speed
    NekDouble a = sqrt(T);

    // Maximum eigenvalue given by |u-a|
    NekDouble C_L = std::max(std::abs(u_L - a), std::abs(u_L + a));
    NekDouble C_R = std::max(std::abs(u_R - a), std::abs(u_R + a));
    NekDouble C = std::max(C_L, C_R);

    // Construct LF flux
    flux[0][j] = 0.5 * (nu_L + nu_R + C * (n_L - n_R));
    flux[1][j] = 0.5 * ((nu_L * u_L + n_L * T) + (nu_R * u_R + n_R * T) +
                        C * (nu_L - nu_R));

    if (is_left) {
      flux[0][j] *= -1.0;
      flux[1][j] *= -1.0;
    }
  }
}
} // namespace NESO::Solvers::H3LAPD

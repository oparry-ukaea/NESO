
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
    : SU::RiemannSolver(session) {
  session->LoadParameter("delta", m_delta, 1.0);
}

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
  boost::ignore_unused(nDim);

  for (int j = 0; j < Fwd[0].size(); ++j) {
    // Field variables
    NekDouble n_L = Fwd[0][j];
    NekDouble nu_L = Fwd[1][j];
    NekDouble n_R = Bwd[0][j];
    NekDouble nu_R = Bwd[1][j];
    NekDouble u_L = nu_L / n_L;

    // Override R conditions for boundaries
    if (j == 0) {
      // x = 0
      n_R = n_L;
      NekDouble ub = u_L + (m_delta - u_L * (-1.0)) * (-1.0);
      nu_R = n_L * (ub);
    } else if (j == 1) {
      n_R = n_L;
      NekDouble ub = u_L + (m_delta - u_L * (1.0)) * (1.0);
      nu_R = n_L * (ub);
    }

    // Velocity
    NekDouble u_R = nu_R / n_R;

    // const temperature
    NekDouble T = 1.0;

    // sound speed
    NekDouble a = m_delta * sqrt(T);

    // Maximum eigenvalue given by |u-a|
    NekDouble C_L = std::max(std::abs(u_L - a), std::abs(u_L + a));
    NekDouble C_R = std::max(std::abs(u_R - a), std::abs(u_R + a));
    NekDouble C = std::max(C_L, C_R);

    // Construct LF flux
    flux[0][j] = 0.5 * (nu_L + nu_R + C * (n_L - n_R));
    flux[1][j] = 0.5 * ((nu_L * u_L + n_L * T) + (nu_R * u_R + n_R * T) +
                        C * (nu_L - nu_R));

    if (j == 0) {
      flux[0][j] *= -1.0;
      flux[1][j] *= -1.0;
    }
  }

  /*
  flux[0][0] = -Fwd[0][0];
  flux[0][1] = Fwd[0][1];
  flux[1][0] = -Fwd[1][0];
  flux[1][1] = Fwd[1][1];
  flux[0][0] = flux[0][1] = flux[1][0] = flux[1][1] = 0.0;
  */
}
} // namespace NESO::Solvers::H3LAPD

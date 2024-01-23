
#include <boost/core/ignore_unused.hpp>

#include "LaxFriedrichsWithOutflowSolver.h"

#include <neso_particles.hpp>

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

  // Load temperature from config
  session->LoadParameter("T", m_T, 1.0);
  // Load mach number from config
  session->LoadParameter("delta", m_mach, 1.0);

  // Compute and store sound speed
  m_cs = std::sqrt(m_T);

  // Compute and store velocity at infinity (outflow velocity)
  // (Sign is accounted for when computing velocity at each boundary)
  m_uinf = m_mach * m_cs;
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

  // // std::cout << "Mesh dim is " << nDim << std::endl;
  // std::cout << "Fwd dims imply " << Fwd.size() << " fields and "
  //           << Fwd[0].size() << " trace points." << std::endl;

  ASSERTL1(CheckScalars("norms"), "Vn not defined.");
  const Nek1DArr &norms = m_scalars["norms"]();

  ASSERTL1(CheckScalars("z"), "z not defined.");
  const Nek1DArr &z = m_scalars["z"]();

  int Nneg = 0, Npos = 0;
  int max_neg = -1, max_pos = -1;
  int min_neg = norms.size(), min_pos = norms.size();
  int z0_nfalse_neg = 0;
  int z2_nfalse_neg = 0;
  for (int j = 0; j < Fwd[0].size(); ++j) {
    // Densities, momenta, velocities in left/fwd and right/bwd trace space
    NekDouble n_L = Fwd[0][j];
    NekDouble nu_L = Fwd[1][j];
    NekDouble n_R = Bwd[0][j];
    NekDouble nu_R = Bwd[1][j];

    if (norms[j] > 0.5) {
      Npos++;
      max_pos = std::max(max_pos, j);
      min_pos = std::min(min_pos, j);
      // std::cout << "norm is " << norms[j] << " at j=" << j << "; z=" << z[j]
      //           << std::endl;
    } else if (norms[j] < -0.5) {
      Nneg++;
      max_neg = std::max(max_neg, j);
      min_neg = std::min(min_neg, j);
    }

    bool is_left = false;
    bool is_right = false;
    switch (nDim) {
    case 1:
      is_left = j == 0;
      is_right = j == 1;
      break;
    case 2:
      NESOASSERT(false, "2D outflow not setup.");
    case 3:
      // Skip transverse boundary (assumes it's also Dirichlet)
      const int nskip = 0;
      const int ntrace_outflow = 10 * 10 * 16;

      if (z[j] < 0.0001) {
        is_left = true;
        if (z[j] > 1e-6 || norms[j] > -0.9999) {
          std::cout << j << ": z=" << z[j] << " with norm " << norms[j]
                    << " was assigned to LEFT boundary!?" << std::endl;
        }
      } else if (z[j] > 1.9999) {
        is_right = true;
        if (z[j] < 1.9999 || norms[j] < 0.9999) {
          std::cout << j << ": z=" << z[j] << " with norm " << norms[j]
                    << " was assigned to LEFT boundary!?" << std::endl;
        }
      } else {
        // if (std::abs(z[j] - 1) > 0.99999 && norms[j] > 0.9999) {
        //   std::cout << j << ": z=" << z[j] << " with norm " << norms[j]
        //             << " is not on boundary!?" << std::endl;
        // }
        if (z[j] < 0.0001) {
          z0_nfalse_neg++;
        } else if (z[j] > 1.9999) {
          z2_nfalse_neg++;
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
    // 0.5*((delta+1/delta)+sqrt((delta+1/delta)^2-4*(z-1)*(z-1)))

    // Velocities in left/fwd and right/bwd trace space
    NekDouble u_L = nu_L / n_L;
    NekDouble u_R = nu_R / n_R;

    // Maximum eigenvalue given by |u-m_cs|; m_cs is the sound speed
    NekDouble C_L = std::max(std::abs(u_L - m_cs), std::abs(u_L + m_cs));
    NekDouble C_R = std::max(std::abs(u_R - m_cs), std::abs(u_R + m_cs));
    NekDouble C = std::max(C_L, C_R);

    // Construct LF flux
    flux[0][j] = 0.5 * (nu_L + nu_R + C * (n_L - n_R));
    flux[1][j] = 0.5 * ((nu_L * u_L + n_L * m_T) + (nu_R * u_R + n_R * m_T) +
                        C * (nu_L - nu_R));

    // Negate left boundary fluxes to account for normal direction
    if (is_left) {
      flux[0][j] *= -1.0;
      flux[1][j] *= -1.0;
    }
  }
  std::cout << z0_nfalse_neg
            << " trace points not identified on the left boundary, "
            << z2_nfalse_neg << " on the right" << std::endl;
  // std::cout << Nneg << " -ve norms between " << min_neg << " and " << max_neg
  //           << std::endl;
  // std::cout << Npos << " +ve norms between " << min_pos << " and " << max_pos
  //           << std::endl
  //           << std::endl;
}
} // namespace NESO::Solvers::H3LAPD

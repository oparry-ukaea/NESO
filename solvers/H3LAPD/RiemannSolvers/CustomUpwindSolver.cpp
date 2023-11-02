
#include <boost/core/ignore_unused.hpp>

#include "CustomUpwindSolver.h"

namespace NESO::Solvers::H3LAPD {

std::string CustomUpwindSolver::solver_name =
    SU::GetRiemannSolverFactory().RegisterCreatorFunction(
        "Custom", CustomUpwindSolver::create, "Custom Upwind solver");

/**
 * @class CustomUpwindSolver
 *
 * @brief Custom upwind scheme Riemann solver.
 *
 * The upwind solver determines the flux based upon an advection
 * velocity \f$\mathbf{V}\f$ and trace normal \f$\mathbf{n}\f$. In
 * particular, the flux for each component of the velocity field is
 * determined by:
 *
 * \f[ \mathbf{f}(u^+,u^-) = \begin{cases} \mathbf{V}u^+, &
 * \mathbf{V}\cdot\mathbf{n} \geq 0,\\ \mathbf{V}u^-, &
 * \mathbf{V}\cdot\mathbf{n} < 0.\end{cases} \f]
 *
 * Here the superscript + and - denotes forwards and backwards spaces
 * respectively.
 */

/**
 * @brief Default constructor.
 */
CustomUpwindSolver::CustomUpwindSolver(
    const LU::SessionReaderSharedPtr &session)
    : SU::RiemannSolver(session) {}

/**
 * @brief Implementation of the upwind solver.
 *
 * The upwind solver assumes that a scalar field Vn is defined, which
 * corresponds with the dot product \f$\mathbf{V}\cdot\mathbf{n}\f$,
 * where \f$\mathbf{V}\f$ is the advection velocity and \f$\mathbf{n}\f$
 * defines the normal of a vertex, edge or face at each quadrature point
 * of the trace space.
 *
 * @param Fwd   Forwards trace space.
 * @param Bwd   Backwards trace space.
 * @param flux  Resulting flux.
 */
void CustomUpwindSolver::v_Solve(const int nDim, const Nek2DArrConstInner &Fwd,
                                 const Nek2DArrConstInner &Bwd,
                                 Nek2DArr &flux) {
  boost::ignore_unused(nDim);

  ASSERTL1(CheckScalars("Foo"), "Foo not defined.");
  const Nek1DArr &foo = m_scalars["Foo"]();

  for (int j = 0; j < foo.size(); ++j) {
    const Nek2DArrConstInner &tmp = foo[j] >= 0 ? Fwd : Bwd;
    // for (int i = 0; i < Fwd.size(); ++i) {
    //   flux[i][j] = foo[j] * tmp[i][j];
    // }
  }
}
} // namespace NESO::Solvers::H3LAPD

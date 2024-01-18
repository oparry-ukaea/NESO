#ifndef H3LAPD_LAXFRIEDRICHSWITHOUTFLOW_H
#define H3LAPD_LAXFRIEDRICHSWITHOUTFLOW_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;

namespace NESO::Solvers::H3LAPD {

typedef Nektar::Array<Nektar::OneD, Nektar::NekDouble> Nek1DArr;
typedef Nektar::Array<Nektar::OneD, Nek1DArr> Nek2DArr;
typedef Nektar::Array<Nektar::OneD, const Nek1DArr> Nek2DArrConstInner;

using Nektar::NekDouble;


/**
 * @brief Subclass Nektar's Riemann solver in order to set up non-standard
 * fluxes.
 */
class LaxFriedrichsWithOutflowSolver : public SU::RiemannSolver {
public:
  SOLVER_UTILS_EXPORT static SU::RiemannSolverSharedPtr
  create(const LU::SessionReaderSharedPtr &pSession) {
    return SU::RiemannSolverSharedPtr(
        new LaxFriedrichsWithOutflowSolver(pSession));
  }

  static std::string solver_name;

protected:
  LaxFriedrichsWithOutflowSolver(const LU::SessionReaderSharedPtr &pSession);

  virtual void v_Solve(const int nDim, const Nek2DArrConstInner &Fwd,
                       const Nek2DArrConstInner &Bwd,
                       Nek2DArr &flux) override final;
};
} // namespace NESO::Solvers::H3LAPD

#endif

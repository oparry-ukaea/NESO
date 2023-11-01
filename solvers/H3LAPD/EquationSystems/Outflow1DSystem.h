#ifndef H3LAPD_1DOUTFLOW_SYSTEM_H
#define H3LAPD_1DOUTFLOW_SYSTEM_H

#include "nektar_interface/utilities.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "DriftReducedSystem.hpp"

namespace NESO::Solvers::H3LAPD {

class Outflow1DSystem : virtual public DriftReducedSystem {
public:
  friend class MemoryManager<Outflow1DSystem>;

  /// Creates an instance of this class.
  static SolverUtils::EquationSystemSharedPtr
  create(const LibUtilities::SessionReaderSharedPtr &session,
         const SpatialDomains::MeshGraphSharedPtr &graph) {
    SolverUtils::EquationSystemSharedPtr p =
        MemoryManager<Outflow1DSystem>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
  }

  /// Name of class.
  static std::string className;

protected:
  Outflow1DSystem(const LibUtilities::SessionReaderSharedPtr &session,
                  const SpatialDomains::MeshGraphSharedPtr &graph);

  virtual void
  explicit_time_int(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, Array<OneD, NekDouble>> &out_arr,
                    const NekDouble time) override;
  virtual void
  get_phi_solve_rhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, NekDouble> &rhs) override;

  virtual void load_params() override;

  void v_InitObject(bool declare_field) override;

private:
  NekDouble m_c;
  NekDouble m_nstar;
  NekDouble m_T;
};

} // namespace NESO::Solvers::H3LAPD
#endif

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

  Array<OneD, NekDouble> m_ztmp;
  Array<OneD, NekDouble> &get_z();

private:
  NekDouble m_c;
  NekDouble m_nstar;
  NekDouble m_T;

  SU::AdvectionSharedPtr m_adv_with_gradP;
  SU::RiemannSolverSharedPtr m_riemann_gradP;

  void
  get_gradP_bulk_flux(const Array<OneD, Array<OneD, NekDouble>> &field_vals,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

  Array<OneD, NekDouble> &get_adv_vel_norm_gradP();
  Array<OneD, NekDouble> &get_trace_normal_z();
};

} // namespace NESO::Solvers::H3LAPD
#endif

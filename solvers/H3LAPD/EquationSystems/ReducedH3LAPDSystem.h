#ifndef H3LAPD_REDUCED_SYSTEM_H
#define H3LAPD_REDUCED_SYSTEM_H

#include "nektar_interface/utilities.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "DriftReducedSystem.hpp"

namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;
namespace NESO::Solvers::H3LAPD {

/**
 * @brief Reduced version of the Hermes-3 LAPD equation system
 */
class ReducedH3LAPDSystem : virtual public DriftReducedSystem {
public:
  friend class MemoryManager<ReducedH3LAPDSystem>;

  /// Creates an instance of this class.
  static SU::EquationSystemSharedPtr
  create(const LU::SessionReaderSharedPtr &session,
         const SD::MeshGraphSharedPtr &graph) {
    SU::EquationSystemSharedPtr p =
        MemoryManager<ReducedH3LAPDSystem>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
  }

  /// Name of class
  static std::string class_name;

protected:
  ReducedH3LAPDSystem(const LU::SessionReaderSharedPtr &session,
                      const SD::MeshGraphSharedPtr &graph);

  void
  add_div_vpar_term(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, Array<OneD, NekDouble>> &out_arr);

  void calc_E_and_adv_vels(
      const Array<OneD, const Array<OneD, NekDouble>> &in_arr) override;

  void
  explicit_time_int(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, Array<OneD, NekDouble>> &out_arr,
                    const NekDouble time) override;

  void
  get_phi_solve_rhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, NekDouble> &rhs) override;

  void v_InitObject(bool declare_field) override;
};

} // namespace NESO::Solvers::H3LAPD
#endif

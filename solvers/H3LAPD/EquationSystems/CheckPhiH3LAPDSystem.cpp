#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "CheckPhiH3LAPDSystem.h"
#include <filesystem>

namespace NESO::Solvers::H3LAPD {
std::string CheckPhiH3LAPDSystem::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "checkPhiLAPD", CheckPhiH3LAPDSystem::create,
        "checkPhi Hermes-3 LAPD equation system");

CheckPhiH3LAPDSystem::CheckPhiH3LAPDSystem(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      DriftReducedSystem(session, graph) {
  // m_required_flds = {"ne", "w", "phi"};
  m_required_flds = {"ne", "Ge", "Gd", "w", "phi", "phi_sln", "phi_diff"};
  m_int_fld_names = {"ne", "w"};
}

void CheckPhiH3LAPDSystem::explicit_time_int(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {

  zero_out_array(out_arr);

  // Solver for electrostatic potential.
  solve_phi(in_arr);

  // Write phi to a separate vtu
  std::filesystem::path phi_output_fname = std::filesystem::path("phi.vtu");
  int phi_idx = m_field_to_index.get_idx("phi");
  if (!std::filesystem::exists(phi_output_fname)) {
    NESO::write_vtu(m_fields[phi_idx], phi_output_fname.string(), "phi");
  }

  // Write phi diff to a separate vtu
  std::filesystem::path phi_diff_output_fname =
      std::filesystem::path("phi_diff.vtu");
  int npts = GetNpoints();
  int phi_diff_idx = m_field_to_index.get_idx("phi_diff");
  int phi_sln_idx = m_field_to_index.get_idx("phi_sln");
  if (!std::filesystem::exists(phi_diff_output_fname)) {
    Vmath::Vsub(npts, m_fields[phi_sln_idx]->GetPhys(), 1,
                m_fields[phi_idx]->GetPhys(), 1,
                m_fields[phi_diff_idx]->UpdatePhys(), 1);
    m_fields[phi_diff_idx]->FwdTrans(m_fields[phi_diff_idx]->GetPhys(),
                                     m_fields[phi_diff_idx]->UpdateCoeffs());
    NESO::write_vtu(m_fields[phi_diff_idx], phi_diff_output_fname.string(),
                    "phi_diff");
  }
}

/**
 * @brief Choose phi solve RHS = w
 *
 * @param in_arr physical values of all fields
 * @param[out] rhs RHS array to pass to Helmsolve
 */
void CheckPhiH3LAPDSystem::get_phi_solve_rhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, NekDouble> &rhs) {
  int npts = GetNpoints();
  int w_idx = m_field_to_index.get_idx("w");
  Vmath::Vcopy(npts, in_arr[w_idx], 1, rhs, 1);
}

/**
 * @brief Post-construction class-initialisation.
 */
void CheckPhiH3LAPDSystem::v_InitObject(bool declare_field) {
  DriftReducedSystem::v_InitObject(declare_field);

  // Bind RHS function for time integration object
  m_ode.DefineOdeRhs(&CheckPhiH3LAPDSystem::explicit_time_int, this);
}

} // namespace NESO::Solvers::H3LAPD

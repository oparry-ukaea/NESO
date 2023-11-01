#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "Outflow1DSystem.h"

namespace NESO::Solvers::H3LAPD {
std::string Outflow1DSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "Outflow1D", Outflow1DSystem::create, "Simple 1D outflow problem");

Outflow1DSystem::Outflow1DSystem(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      DriftReducedSystem(session, graph) {
  m_required_flds = {"ne", "Ge"};
  m_int_fld_names = {"ne", "Ge"};
}

void Outflow1DSystem::explicit_time_int(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {

  // Check inarray for NaNs
  for (auto &var : {"ne", "Ge"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < inarray[fidx].size(); ii++) {
      if (!std::isfinite(inarray[fidx][ii])) {
        std::cout << "NaN in field " << var << ", aborting." << std::endl;
        exit(1);
      }
    }
  }

  zero_out_array(outarray);

  // Compute ve from Ge

  // Get field indices
  int nPts = GetNpoints();
  int ne_idx = m_field_to_index.get_idx("ne");
  int Ge_idx = m_field_to_index.get_idx("Ge");

  // v_par,e = Ge / max(ne,n_floor) / me
  for (auto idim = 0; idim < m_graph->GetSpaceDimension(); idim++) {
    for (auto ii = 0; ii < nPts; ii++) {
      m_adv_vel_elec[idim][ii] =
          inarray[Ge_idx][ii] /
          std::max(inarray[ne_idx][ii], m_n_ref * m_n_floor_fac);
    }
  }

  add_adv_terms({"ne"}, m_adv_elec, m_ExB_vel, inarray, outarray, time);

  // // Advect both ne and w using ExB velocity
  // add_adv_terms({"ne"}, m_advElec, m_vExB, inarray, outarray, time);
  // add_adv_terms({"w"}, m_advVort, m_vExB, inarray, outarray, time);

  // // Add \alpha*(\phi-n_e) to RHS
  // Array<OneD, NekDouble> HWterm_2D_alpha(nPts);
  // Vmath::Vsub(nPts, m_fields[phi_idx]->GetPhys(), 1,
  //             m_fields[ne_idx]->GetPhys(), 1, HWterm_2D_alpha, 1);
  // Vmath::Smul(nPts, m_alpha, HWterm_2D_alpha, 1, HWterm_2D_alpha, 1);
  // Vmath::Vadd(nPts, outarray[w_idx], 1, HWterm_2D_alpha, 1, outarray[w_idx],
  // 1); Vmath::Vadd(nPts, outarray[ne_idx], 1, HWterm_2D_alpha, 1,
  // outarray[ne_idx],
  //             1);

  // // Add \kappa*\dpartial\phi/\dpartial y to RHS
  // Array<OneD, NekDouble> HWterm_2D_kappa(nPts);
  // m_fields[phi_idx]->PhysDeriv(1, m_fields[phi_idx]->GetPhys(),
  //                              HWterm_2D_kappa);
  // Vmath::Vsub(nPts, outarray[ne_idx], 1, HWterm_2D_kappa, 1,
  // outarray[ne_idx],
  //             1);
}

void Outflow1DSystem::load_params() {
  DriftReducedSystem::load_params();

  // c
  m_session->LoadParameter("c", m_c);

  // nstar
  m_session->LoadParameter("n_star", m_nstar);

  // T
  m_session->LoadParameter("T", m_T);
}

// dummy
void Outflow1DSystem::get_phi_solve_rhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, NekDouble> &rhs) {
  // Do nothing
}

/**
 * @brief Post-construction class-initialisation.
 */
void Outflow1DSystem::v_InitObject(bool declare_field) {
  DriftReducedSystem::v_InitObject(declare_field);

  // Bind RHS function for time integration object
  m_ode.DefineOdeRhs(&Outflow1DSystem::explicit_time_int, this);
}

} // namespace NESO::Solvers::H3LAPD

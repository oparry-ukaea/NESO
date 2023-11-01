#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "ReducedH3LAPDSystem.h"

namespace NESO::Solvers::H3LAPD {
std::string ReducedH3LAPDSystem::class_name =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "reducedLAPD", ReducedH3LAPDSystem::create,
        "Reduced Hermes-3 LAPD equation system");

ReducedH3LAPDSystem::ReducedH3LAPDSystem(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      DriftReducedSystem(session, graph) {
  m_required_flds = {"ne", "w", "phi"};
}

void ReducedH3LAPDSystem::add_div_vpar_term(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr) {
  int ne_idx = m_field_to_index.get_idx("ne");
  int w_idx = m_field_to_index.get_idx("w");
  int npts = GetNpoints();
  Array<OneD, NekDouble> div_nvpar_term(npts);
  Vmath::Vmul(npts, in_arr[ne_idx], 1, m_par_vel_elec, 1, div_nvpar_term, 1);
  m_fields[ne_idx]->PhysDeriv(2, div_nvpar_term, div_nvpar_term);
  Vmath::Vsub(npts, out_arr[w_idx], 1, div_nvpar_term, 1, out_arr[w_idx], 1);
}

/**
 * @brief Compute E = \f$ -\nabla\phi\f$, \f$ v_{E\times B}\f$ and the advection
 * velocities used in the ne/Ge, Gd equations.
 * @param in_arr array of field physvals
 */
void ReducedH3LAPDSystem::calc_E_and_adv_vels(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr) {
  DriftReducedSystem::calc_E_and_adv_vels(in_arr);
  // int phi_idx = m_field_to_index.get_idx("phi");
  int npts = GetNpoints();

  // vpar
  int ne_idx = m_field_to_index.get_idx("ne");
  Array<OneD, NekDouble> tmpx(npts), tmpy(npts), tmpz(npts);
  m_fields[ne_idx]->GetCoords(tmpx, tmpy, tmpz);
  m_session->GetFunction("vpar", ne_idx)
      ->Evaluate(tmpx, tmpy, tmpz, m_par_vel_elec);
  for (auto iDim = 0; iDim < m_graph->GetSpaceDimension(); iDim++) {
    Vmath::Svtvp(npts, m_b_unit[iDim], m_par_vel_elec, 1, m_ExB_vel[iDim], 1,
                 m_adv_vel_elec[iDim], 1);
  }
}

void ReducedH3LAPDSystem::explicit_time_int(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {

  for (auto &var : {"ne", "Gd", "Ge", "w"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < in_arr[fidx].size(); ii++) {
      std::stringstream err_msg;
      err_msg << "Found NaN in field " << var;
      NESOASSERT(std::isfinite(in_arr[fidx][ii]), err_msg.str().c_str());
    }
  }

  zero_out_array(out_arr);

  // Solver for electrostatic potential.
  solve_phi(in_arr);

  // Calculate electric field from Phi, and compute v_ExB
  calc_E_and_adv_vels(in_arr);

  // Add ne advection term to out_arr
  add_adv_terms({"ne"}, m_adv_elec, m_adv_vel_elec, in_arr, out_arr, time);
  // Add w advection term to out_arr
  add_adv_terms({"w"}, m_adv_vort, m_ExB_vel, in_arr, out_arr, time);

  add_div_vpar_term(in_arr, out_arr);
}

/**
 * @brief Choose phi solve RHS = w
 *
 * @param in_arr physical values of all fields
 * @param[out] rhs RHS array to pass to Helmsolve
 */
void ReducedH3LAPDSystem::get_phi_solve_rhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, NekDouble> &rhs) {
  int npts = GetNpoints();
  int w_idx = m_field_to_index.get_idx("w");
  Vmath::Vcopy(npts, in_arr[w_idx], 1, rhs, 1);
}

/**
 * @brief Post-construction class-initialisation.
 */
void ReducedH3LAPDSystem::v_InitObject(bool declare_field) {
  DriftReducedSystem::v_InitObject(declare_field);

  // Bind RHS function for time integration object
  m_ode.DefineOdeRhs(&ReducedH3LAPDSystem::explicit_time_int, this);
}

} // namespace NESO::Solvers::H3LAPD

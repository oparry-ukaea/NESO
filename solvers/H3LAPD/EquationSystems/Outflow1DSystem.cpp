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
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {

  // Check in_arr for NaNs
  for (auto &var : {"ne", "Ge"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < in_arr[fidx].size(); ii++) {
      std::stringstream tmp;
      tmp << "NaN in field " << var << ", aborting." << std::endl;
      NESOASSERT(std::isfinite(in_arr[fidx][ii]), tmp.str().c_str());
    }
  }

  zero_out_array(out_arr);

  // Get field indices
  int nPts = GetNpoints();
  int ne_idx = m_field_to_index.get_idx("ne");
  int Ge_idx = m_field_to_index.get_idx("Ge");

  // Set advection velocity from parallel momentum and density, accounting for
  // density floor
  for (auto idim = 0; idim < m_graph->GetSpaceDimension(); idim++) {
    for (auto ii = 0; ii < nPts; ii++) {
      m_adv_vel_elec[idim][ii] =
          in_arr[Ge_idx][ii] /
          std::max(in_arr[ne_idx][ii], m_n_ref * m_n_floor_fac);
    }
  }

  // Add advection terms to ne and Ge equations, gradP flux to Ge equation
  add_adv_terms({"ne", "Ge"}, m_adv_with_gradP, m_adv_vel_elec, in_arr, out_arr,
                time);
  // add_adv_terms({"ne", "Ge"}, m_adv_elec, m_adv_vel_elec, in_arr, out_arr,
  //               time);

  // Add source to ne equation
  add_density_source(out_arr);
}

void Outflow1DSystem::load_params() {
  DriftReducedSystem::load_params();

  // No additional params to load
}

/**
 * No phi solve for this system - dummy function to keep the compiler happy.
 */
void Outflow1DSystem::get_phi_solve_rhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, NekDouble> &rhs) {
  // Do nothing
}

void Outflow1DSystem::get_gradP_bulk_flux(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  NESOASSERT(flux.size() == 2, "Expecting flux array with outer dim size 2");
  // Add standard bulk flux terms
  get_flux_vector(field_vals, m_adv_vel_elec, flux);

  /* Code to add gradP contribution assumes fields are ne,Ge.
  Make sure that's the case
   */
  int ne_idx = m_field_to_index.get_idx("ne");
  int Ge_idx = m_field_to_index.get_idx("Ge");
  int flow_dim = m_graph->GetSpaceDimension() - 1;
  NESOASSERT(ne_idx == 0 && Ge_idx == 1, "Unexpected indices");
  NESOASSERT(flux[Ge_idx].size() == m_graph->GetSpaceDimension(),
             "Unexpected flux dim");

  // Add gradP term
  for (int i = 0; i < field_vals[ne_idx].size(); ++i) {
    flux[Ge_idx][flow_dim][i] += field_vals[ne_idx][i];
  }
}

/**
 * @brief Compute trace-normal advection velocities
 */
Array<OneD, NekDouble> &Outflow1DSystem::get_adv_vel_norm_gradP() {
  Array<OneD, NekDouble> norm_vels(GetTraceNpoints());
  return get_adv_vel_norm(norm_vels, m_adv_vel_elec);
}

Array<OneD, NekDouble> &Outflow1DSystem::get_trace_normal_z() {
  if (m_graph->GetSpaceDimension() == 3) {
    return m_traceNormals[2];
  } else if (m_graph->GetSpaceDimension() == 1) {
    return m_traceNormals[0];
  } else {
    NESOASSERT(false, "Not set up for 2D");
    // Keep compiler happy
    return m_traceNormals[0];
  }
}

/**
 * @brief Post-construction class-initialisation.
 */
void Outflow1DSystem::v_InitObject(bool declare_field) {
  DriftReducedSystem::v_InitObject(declare_field);

  // Bind RHS function for time integration object
  m_ode.DefineOdeRhs(&Outflow1DSystem::explicit_time_int, this);

  // Setup custom riemann solver and quasi-advection object for pressure
  // gradient term
  //

  m_adv_with_gradP =
      SU::GetAdvectionFactory().CreateInstance(m_adv_type, m_adv_type);
  // Set callback function to compute bulk flux
  m_adv_with_gradP->SetFluxVector(&Outflow1DSystem::get_gradP_bulk_flux, this);

  // Create Riemann solver
  m_riemann_gradP =
      SU::GetRiemannSolverFactory().CreateInstance("Custom", m_session);
  m_riemann_gradP->SetScalar("nz", &Outflow1DSystem::get_trace_normal_z, this);
  m_riemann_gradP->SetScalar("Vn", &Outflow1DSystem::get_adv_vel_norm_gradP,
                             this);

  // Bind to quasi-advection object
  m_adv_with_gradP->SetRiemannSolver(m_riemann_gradP);
  m_adv_with_gradP->InitObject(m_session, m_fields);
}

} // namespace NESO::Solvers::H3LAPD

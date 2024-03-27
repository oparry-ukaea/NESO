#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include <cmath>
#include <mpi.h>
#include <random>

#include "DriftReducedSystem.hpp"

namespace NESO::Solvers::H3LAPD {
DriftReducedSystem::DriftReducedSystem(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), AdvectionSystem(session, graph),
      m_field_to_index(session->GetVariables()),
      m_adv_vel_elec(graph->GetSpaceDimension()),
      m_ExB_vel(graph->GetSpaceDimension()), m_E(graph->GetSpaceDimension()) {
  // Construct particle system
  m_particle_sys = std::make_shared<NeutralParticleSystem>(session, graph);
}

/**
 * @brief Compute advection terms and add them to an output array
 * @details For each field listed in @p field_names copy field pointers from
 * m_fields and physical values from @p in_arr to temporary arrays and create a
 * temporary output array of the same size. Call Advect() on @p adv_obj passing
 * the temporary arrays. Finally, loop over @p eqn_labels to determine which
 * element(s) of @p out_arr to subtract the results from.
 *
 * N.B. The creation of temporary arrays is necessary to bypass restrictions in
 * the Nektar advection API.
 *
 * @param field_names List of field names to compute advection terms for
 * @param adv_obj Nektar advection object
 * @param adv_vel Array of advection velocities (outer dim size = nfields)
 * @param in_arr Physical values for *all* fields
 * @param[out] out_arr RHS array (for *all* fields)
 * @param time Simulation time
 * @param eqn_labels List of field names identifying indices in out_arr to add
 * the result to. Defaults to @p field_names
 */
void DriftReducedSystem::add_adv_terms(
    std::vector<std::string> field_names, const SU::AdvectionSharedPtr adv_obj,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time,
    std::vector<std::string> eqn_labels) {

  // Default is to add result of advecting field f to the RHS of df/dt equation
  if (eqn_labels.empty()) {
    eqn_labels = std::vector(field_names);
  } else {
    ASSERTL1(field_names.size() == eqn_labels.size(),
             "add_adv_terms: Number of quantities being advected must match "
             "the number of equation labels.");
  }

  int nfields = field_names.size();
  int npts = GetNpoints();

  /* Make temporary copies of target fields, in_arr vals and initialise a
   * temporary output array
   */
  Array<OneD, MR::ExpListSharedPtr> tmp_fields(nfields);
  Array<OneD, Array<OneD, NekDouble>> tmp_inarray(nfields);
  Array<OneD, Array<OneD, NekDouble>> tmp_outarray(nfields);
  for (auto ii = 0; ii < nfields; ii++) {
    int idx = m_field_to_index.get_idx(field_names[ii]);
    tmp_fields[ii] = m_fields[idx];
    tmp_inarray[ii] = Array<OneD, NekDouble>(npts);
    Vmath::Vcopy(npts, in_arr[idx], 1, tmp_inarray[ii], 1);
    tmp_outarray[ii] = Array<OneD, NekDouble>(out_arr[idx].size());
  }
  // Compute advection terms; result is returned in temporary output array
  adv_obj->Advect(tmp_fields.size(), tmp_fields, adv_vel, tmp_inarray,
                  tmp_outarray, time);

  // Subtract temporary output array from the appropriate indices of out_arr
  for (auto ii = 0; ii < nfields; ii++) {
    int idx = m_field_to_index.get_idx(eqn_labels[ii]);
    Vmath::Vsub(out_arr[idx].size(), out_arr[idx], 1, tmp_outarray[ii], 1,
                out_arr[idx], 1);
  }
}

/**
 * @brief Add (density) source term via a Nektar session function.
 *
 * @details Looks for a function called "dens_src", evaluates it, and adds the
 * result to @p out_arr
 *
 * @param[out] out_arr RHS array to add the source too
 * @todo Check function exists, rather than relying on Nektar ASSERT
 */
void DriftReducedSystem::add_density_source(
    Array<OneD, Array<OneD, NekDouble>> &out_arr) {

  int ne_idx = m_field_to_index.get_idx("ne");
  int npts = GetNpoints();
  Array<OneD, NekDouble> tmpx(npts), tmpy(npts), tmpz(npts);
  m_fields[ne_idx]->GetCoords(tmpx, tmpy, tmpz);
  Array<OneD, NekDouble> dens_src(npts, 0.0);
  LU::EquationSharedPtr dens_src_func =
      m_session->GetFunction("dens_src", ne_idx);
  dens_src_func->Evaluate(tmpx, tmpy, tmpz, dens_src);
  Vmath::Vadd(npts, out_arr[ne_idx], 1, dens_src, 1, out_arr[ne_idx], 1);
}

/**
 * @brief Adds particle sources.
 * @details For each <field_name> in @p target_fields , look for another field
 * called <field_name>_src. If it exists, add the physical values of
 * field_name_src to the appropriate index of @p out_arr.
 *
 * @param target_fields list of Nektar field names for which to look for a
 * '_src' counterpart
 *  @param[out] out_arr      the RHS array
 *
 */
void DriftReducedSystem::add_particle_sources(
    std::vector<std::string> target_fields,
    Array<OneD, Array<OneD, NekDouble>> &out_arr) {
  for (auto target_field : target_fields) {
    int src_field_idx = m_field_to_index.get_idx(target_field + "_src");

    if (src_field_idx >= 0) {
      // Check that the target field is one that is time integrated
      auto tmp_it = std::find(m_int_fld_names.cbegin(), m_int_fld_names.cend(),
                              target_field);
      ASSERTL0(tmp_it != m_int_fld_names.cend(),
               "Target for particle source ['" + target_field +
                   "'] does not appear in the list of time-integrated fields "
                   "(m_int_fld_names).")
      /*
      N.B. out_arr can be smaller than m_fields if any fields aren't
      time-integrated, so can't just use out_arr_idx =
      m_field_to_index.get_idx(target_field)
       */
      auto out_arr_idx = std::distance(m_int_fld_names.cbegin(), tmp_it);
      Vmath::Vadd(out_arr[out_arr_idx].size(), out_arr[out_arr_idx], 1,
                  m_fields[src_field_idx]->GetPhys(), 1, out_arr[out_arr_idx],
                  1);
    }
  }
}

/**
 *  @brief Compute E = \f$ -\nabla\phi\f$, \f$ v_{E\times B}\f$ and the
 * advection velocities used in the ne/Ge, Gd equations.
 *
 * @param in_arr array of field phys vals
 */
void DriftReducedSystem::calc_E_and_adv_vels(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr) {
  int phi_idx = m_field_to_index.get_idx("phi");
  int npts = GetNpoints();
  m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(), m_E[0], m_E[1],
                               m_E[2]);
  Vmath::Neg(npts, m_E[0], 1);
  Vmath::Neg(npts, m_E[1], 1);
  Vmath::Neg(npts, m_E[2], 1);

  // v_ExB = Evec x Bvec / B^2
  Vmath::Svtsvtp(npts, m_B[2] / m_Bmag / m_Bmag, m_E[1], 1,
                 -m_B[1] / m_Bmag / m_Bmag, m_E[2], 1, m_ExB_vel[0], 1);
  Vmath::Svtsvtp(npts, m_B[0] / m_Bmag / m_Bmag, m_E[2], 1,
                 -m_B[2] / m_Bmag / m_Bmag, m_E[0], 1, m_ExB_vel[1], 1);
  Vmath::Svtsvtp(npts, m_B[1] / m_Bmag / m_Bmag, m_E[0], 1,
                 -m_B[0] / m_Bmag / m_Bmag, m_E[1], 1, m_ExB_vel[2], 1);
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * @details This routine projects the @p in_arr input and ensures the @p
 * out_arr output lives in the correct space. Since we are hard-coding DG, this
 * corresponds to a simple copy from in to out, since no elemental
 * connectivity is required and the output of the RHS function is
 * polynomial.
 *
 * @param in_arr Unprojected values
 * @param[out] out_arr Projected values
 * @param time Current simulation time
 *
 */
void DriftReducedSystem::do_ode_projection(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {
  int num_vars = in_arr.size();
  int npoints = in_arr[0].size();

  for (int i = 0; i < num_vars; ++i) {
    Vmath::Vcopy(npoints, in_arr[i], 1, out_arr[i], 1);
  }
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel        Advection velocities for each field
 */
Array<OneD, NekDouble> &DriftReducedSystem::get_adv_vel_norm(
    Array<OneD, NekDouble> &trace_vel_norm,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel) {
  // Number of trace (interface) points
  int num_trace_pts = GetTraceNpoints();
  // Auxiliary variable to compute normal velocities
  Array<OneD, NekDouble> tmp(num_trace_pts);

  // Zero previous values
  Vmath::Zero(num_trace_pts, trace_vel_norm, 1);

  //  Compute dot product of advection velocity with the trace normals and store
  for (int i = 0; i < adv_vel.size(); ++i) {
    m_fields[0]->ExtractTracePhys(adv_vel[i], tmp);
    Vmath::Vvtvp(num_trace_pts, m_traceNormals[i], 1, tmp, 1, trace_vel_norm, 1,
                 trace_vel_norm, 1);
  }
  return trace_vel_norm;
}

/**
 *  @brief Compute trace-normal advection velocities for the electron density.
 */
Array<OneD, NekDouble> &DriftReducedSystem::get_adv_vel_norm_elec() {
  return get_adv_vel_norm(m_norm_vel_elec, m_adv_vel_elec);
}

/**
 * @brief Compute trace-normal advection velocities for the vorticity equation.
 */
Array<OneD, NekDouble> &DriftReducedSystem::get_adv_vel_norm_vort() {
  return get_adv_vel_norm(m_norm_vel_vort, m_ExB_vel);
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param  adv_vel    Advection velocities for each advection field
 * @param[out] flux       Flux array
 */
void DriftReducedSystem::get_flux_vector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  ASSERTL1(flux[0].size() == adv_vel.size(),
           "Dimension of flux array and advection velocity array do not match");
  int npts = field_vals[0].size();

  for (auto i = 0; i < flux.size(); ++i) {
    for (auto j = 0; j < flux[0].size(); ++j) {
      Vmath::Vmul(npts, field_vals[i], 1, adv_vel[j], 1, flux[i][j], 1);
    }
  }
}

/**
 * @brief Construct the flux vector for the diffusion problem.
 * @todo not implemented
 */
void DriftReducedSystem::get_flux_vector_diff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor) {
  std::cout << "*** GetFluxVectorDiff not defined! ***" << std::endl;
}

/**
 * @brief Compute the flux vector for advection in the electron density and
 * momentum equations.
 *
 * @param field_vals   Array of Fields ptrs
 * @param[out] flux         Resulting flux array
 */
void DriftReducedSystem::get_flux_vector_elec(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  get_flux_vector(field_vals, m_adv_vel_elec, flux);
}

/**
 * @brief Compute the flux vector for advection in the vorticity equation.
 *
 * @param field_vals   Array of Fields ptrs
 * @param[out] flux        Resulting flux array
 */
void DriftReducedSystem::get_flux_vector_vort(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  // Advection velocity is v_ExB in the vorticity equation
  get_flux_vector(field_vals, m_ExB_vel, flux);
}

/**
 * @brief Load all required session parameters into member variables.
 */
void DriftReducedSystem::load_params() {
  // Type of advection to use -- in theory we also support flux reconstruction
  // for quad-based meshes, or you can use a standard convective term if you
  // were fully continuous in space. Default is DG.
  m_session->LoadSolverInfo("AdvectionType", m_adv_type, "WeakDG");

  // ***Assumes field aligned with z-axis***
  // Magnetic field strength. Fix B = [0, 0, Bxy] for now
  m_B = std::vector<NekDouble>(m_graph->GetSpaceDimension(), 0);
  m_session->LoadParameter("Bxy", m_B[2], 0.1);

  // Coefficient factors for potential solve
  m_session->LoadParameter("d00", m_d00, 1);
  m_session->LoadParameter("d11", m_d11, 1);
  m_session->LoadParameter("d22", m_d22, 1);

  // Factor to set density floor; default to 1e-5 (Hermes-3 default)
  m_session->LoadParameter("n_floor_fac", m_n_floor_fac, 1e-5);

  // Reference number density
  m_session->LoadParameter("nRef", m_n_ref, 1.0);

  // Type of Riemann solver to use. Default = "Upwind"
  m_session->LoadSolverInfo("UpwindType", m_riemann_solver_type, "Upwind");

  // Particle-related parameters
  m_session->LoadParameter("num_particle_steps_per_fluid_step",
                           m_num_part_substeps, 1);
  m_session->LoadParameter("particle_num_write_particle_steps",
                           m_num_write_particle_steps, 0);
  m_part_timestep = m_timestep / m_num_part_substeps;

  // Compute some properties derived from params
  m_Bmag = std::sqrt(m_B[0] * m_B[0] + m_B[1] * m_B[1] + m_B[2] * m_B[2]);
  m_b_unit = std::vector<NekDouble>(m_graph->GetSpaceDimension());
  for (auto idim = 0; idim < m_b_unit.size(); idim++) {
    m_b_unit[idim] = (m_Bmag > 0) ? m_B[idim] / m_Bmag : 0.0;
  }
}

void DriftReducedSystem::get_scaled_coords(
    Array<OneD, Array<OneD, NekDouble>> &coords) {
  int nPts = GetNpoints();
  ASSERTL1(
      coords.size() == 3,
      "get_scaled_coords: supplied coord array doesn't have 3 dimensions.");
  ASSERTL1(coords[0].size() == nPts,
           "get_scaled_coords: supplied coord array doesn't have the right "
           "number of points.");

  // Local mesh extent
  auto vertices = m_graph->GetAllPointGeoms();
  double origin[3];
  double extent[3];
  for (int dimx = 0; dimx < 3; dimx++) {
    origin[dimx] = std::numeric_limits<double>::max();
    extent[dimx] = std::numeric_limits<double>::min();
  }

  for (auto &vx : vertices) {
    Nektar::NekDouble x, y, z;
    vx.second->GetCoords(x, y, z);
    origin[0] = std::min(origin[0], x);
    origin[1] = std::min(origin[1], y);
    origin[2] = std::min(origin[2], z);
    extent[0] = std::max(extent[0], x);
    extent[1] = std::max(extent[1], y);
    extent[2] = std::max(extent[2], z);
  }

  // Global mesh extent
  NekDouble global_origin[3];
  NekDouble global_extent[3];
  MPICHK(MPI_Allreduce(origin, global_origin, 3, MPI_DOUBLE, MPI_MIN,
                       MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(extent, global_extent, 3, MPI_DOUBLE, MPI_MAX,
                       MPI_COMM_WORLD));

  // Scale coords
  // x between 0 and 1; y, z between 0 and 2pi
  for (int dim = 0; dim < 3; dim++) {
    NekDouble scale;
    if (dim == 0) {
      scale = 1.0;
    } else {
      scale = 2 * M_PI;
    }
    global_extent[dim] -= global_origin[dim];
    Vmath::Sadd(nPts, -1 * global_origin[dim], coords[dim], 1, coords[dim], 1);
    Vmath::Smul(nPts, scale / global_extent[dim], coords[dim], 1, coords[dim],
                1);
  }
}

double eval_mode(double x, int mode, double phase) {
  double prefac = 1.0 / (1.0 + std::abs(mode - 4.0));
  return prefac * prefac * std::cos(mode * x + phase);
}

void DriftReducedSystem::get_phases(std::vector<NekDouble> &phases,
                                    int n_modes) {
  auto comm = m_session->GetComm();
  if (m_session->GetComm()->TreatAsRankZero()) {
    int seed = 1;
    // Generator
    std::mt19937 rng(seed);
    // Distribution: uniform on [-pi, +pi]
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);

    for (int imode = 0; imode < n_modes; imode++) {
      phases[imode] = dist(rng);
    }
    comm->Bcast(phases, 0);
  } else {
    comm->Bcast(phases, 0);
  }
}

void DriftReducedSystem::mixmode(const Array<OneD, NekDouble> &xarr,
                                 Array<OneD, NekDouble> &result, int n_modes) {
  // Set phases on root and broadcast to other ranks
  std::vector<NekDouble> phases(n_modes);
  get_phases(phases, n_modes);

  Vmath::Zero(result.size(), result, 1);
  for (int iPt = 0; iPt < xarr.size(); iPt++) {
    for (int imode = 0; imode < n_modes; imode++) {
      result[iPt] += eval_mode(xarr[iPt], imode, phases[imode]);
    }
  }
}

void DriftReducedSystem::set_mixed_mode_ICs(const double &scale,
                                            const int &n_modes,
                                            const int &peak_mode) {
  int nPts = GetNpoints();
  // Get scaled coords
  Array<OneD, Array<OneD, NekDouble>> coords(3);
  for (auto ii = 0; ii < coords.size(); ii++) {
    coords[ii] = Array<OneD, NekDouble>(nPts);
  }
  m_fields[0]->GetCoords(coords[0], coords[1], coords[2]);
  get_scaled_coords(coords);

  // Compute 2*PI*x and y-z
  Array<OneD, NekDouble> twoPIx(nPts), yminusz(nPts);
  Vmath::Smul(nPts, 2 * M_PI, coords[0], 1, twoPIx, 1);
  Vmath::Vsub(nPts, coords[1], 1, coords[2], 1, yminusz, 1);

  std::vector<std::string> mixed_mode_vars = {"w"};
  for (auto &var : mixed_mode_vars) {
    int var_idx = m_field_to_index.get_idx(var);
    // Get non-const ref to var phys vals
    auto phys_vals = m_fields[var_idx]->UpdatePhys();

    // Set phys vals = scale * mixmode(2*pi*x) * mixmode(z - y)
    Array<OneD, NekDouble> term1(nPts), term2(nPts);
    mixmode(twoPIx, term1, n_modes);
    mixmode(yminusz, term2, n_modes);
    Vmath::Vmul(nPts, term1, 1, term2, 1, phys_vals, 1);
    std::cout << "mixmode scale is " << scale << std::endl;
    Vmath::Smul(nPts, scale, phys_vals, 1, phys_vals, 1);
    // Update coeffs
    m_fields[var_idx]->FwdTrans(phys_vals, m_fields[var_idx]->UpdateCoeffs());
  }
}

/**
 * @brief Utility function to print the size of a 1D Nektar array.
 * @param arr Array to print the size of
 * @param label Label to include in the output message
 * @param all_tasks If true, print message on all tasks, else print only on task
 * 0 (default=false)
 */
void DriftReducedSystem::print_arr_size(const Array<OneD, NekDouble> &arr,
                                        std::string label, bool all_tasks) {
  if (m_session->GetComm()->TreatAsRankZero() || all_tasks) {
    if (!label.empty()) {
      std::cout << label << " ";
    }
    std::cout << "size = " << arr.size() << std::endl;
  }
}

/**
 * @brief Utility function to print values in a 1D Nektar array.
 *
 * @param arr Nektar array from which to extract values
 * @param num number of values to report
 * @param stride stride between indices (first value has index 0)
 * @param label label to use for the array when reporting values
 * @param all_tasks flag to output the result on all tasks (default is just task
 * 0)
 */
void DriftReducedSystem::print_arr_vals(const Array<OneD, NekDouble> &arr,
                                        int num, int stride, std::string label,
                                        bool all_tasks) {
  if (m_session->GetComm()->TreatAsRankZero() || all_tasks) {
    if (!label.empty()) {
      std::cout << "[" << label << "]" << std::endl;
    }
    int ii_max = std::min(static_cast<int>(arr.size()), num * stride);
    for (auto ii = 0; ii < ii_max; ii = ii + stride) {
      std::cout << "  " << std::setprecision(12) << arr[ii] << std::endl;
    }
  }
}

/**
 * @brief Calls HelmSolve to solve for the electric potential, given the
 * right-hand-side returned by get_phi_solve_rhs
 *
 * @param in_arr Array of physical field values
 */
void DriftReducedSystem::solve_phi(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr) {

  // Field indices
  int npts = GetNpoints();
  int phi_idx = m_field_to_index.get_idx("phi");

  // Define rhs
  Array<OneD, NekDouble> rhs(npts);
  get_phi_solve_rhs(in_arr, rhs);

  // Set up factors for electrostatic potential solve
  StdRegions::ConstFactorMap factors;
  // Helmholtz => Poisson (lambda = 0)
  factors[StdRegions::eFactorLambda] = 0.0;
  // Set coefficient factors
  factors[StdRegions::eFactorCoeffD00] = m_d00;
  factors[StdRegions::eFactorCoeffD11] = m_d11;
  factors[StdRegions::eFactorCoeffD22] = m_d22;

  // Solve for phi. Output of this routine is in coefficient (spectral)
  // space, so backwards transform to physical space since we'll need that
  // for the advection step & computing drift velocity.
  m_fields[phi_idx]->HelmSolve(rhs, m_fields[phi_idx]->UpdateCoeffs(), factors);
  m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(),
                              m_fields[phi_idx]->UpdatePhys());
}

/**
 * @brief Check required fields are all defined and have the same number of quad
 * points
 */
void DriftReducedSystem::validate_fields() {
  int npts_exp = GetNpoints();
  for (auto &fld_name : m_required_flds) {
    int idx = m_field_to_index.get_idx(fld_name);
    // Check field exists
    ASSERTL0(idx >= 0, "Required field [" + fld_name + "] is not defined.");
    // Check fields all have the same number of quad points
    int npts = m_fields[idx]->GetNpoints();
    ASSERTL0(npts == npts_exp,
             "Expecting " + std::to_string(npts_exp) +
                 " quad points, but field '" + fld_name + "' has " +
                 std::to_string(npts) +
                 ". Check NUMMODES is the same for all required fields.");
  }
}

void DriftReducedSystem::v_SetInitialConditions(NekDouble initialtime,
                                                bool dumpInitialConditions,
                                                const int domain) {

  boost::ignore_unused(initialtime);

  if (m_session->GetComm()->GetRank() == 0) {
    cout << "Initial Conditions:" << endl;
  }

  if (m_session->DefinesFunction("InitialConditions")) {
    GetFunction("InitialConditions")
        ->Evaluate(m_session->GetVariables(), m_fields, m_time, domain);
    // Enforce C0 Continutiy of initial condiiton
    if ((m_projectionType == MultiRegions::eGalerkin) ||
        (m_projectionType == MultiRegions::eMixed_CG_Discontinuous)) {
      for (int i = 0; i < m_fields.size(); ++i) {
        m_fields[i]->LocalToGlobal();
        m_fields[i]->GlobalToLocal();
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                              m_fields[i]->UpdatePhys());
      }
    }

    if (m_session->GetComm()->GetRank() == 0) {

      for (int i = 0; i < m_fields.size(); ++i) {
        std::string varName = m_session->GetVariable(i);
        cout << "  - Field " << varName << ": "
             << GetFunction("InitialConditions")->Describe(varName, domain)
             << endl;
      }
    }
  } else {
    int nq = m_fields[0]->GetNpoints();
    for (int i = 0; i < m_fields.size(); i++) {
      Vmath::Zero(nq, m_fields[i]->UpdatePhys(), 1);
      m_fields[i]->SetPhysState(true);
      Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
      if (m_session->GetComm()->GetRank() == 0) {
        cout << "  - Field " << m_session->GetVariable(i) << ": 0 (default)"
             << endl;
      }
    }
  }
  if (m_session->DefinesParameter("mixmode_ICs")) {
    // Set 'mixed mode' distribution, a la BOUT++
    //   Add <n_modes> in total

    NekDouble scale_fudge;
    m_session->LoadParameter("mixmode_scale_fudge", scale_fudge, 1.0);

    int n_modes;
    m_session->LoadParameter("mixmode_nmodes", n_modes, 14);
    //   Max amplitude is at <peak_mode>
    const int peak_mode = 4;
    const double fluctuation_scale = 0.1 * scale_fudge;
    set_mixed_mode_ICs(fluctuation_scale, n_modes, peak_mode);
  }

  // Compute init phi
  Array<OneD, Array<OneD, NekDouble>> phys_vals(m_fields.size());
  int npts = GetNpoints();
  for (auto ii = 0; ii < phys_vals.size(); ii++) {
    phys_vals[ii] = Array<OneD, NekDouble>(npts);
    Vmath::Vcopy(npts, m_fields[ii]->GetPhys(), 1, phys_vals[ii], 1);
  }
  solve_phi(phys_vals);

  if (dumpInitialConditions && m_checksteps && m_nchk == 0) {
    Checkpoint_Output(m_nchk);
  }
  ++m_nchk;
}

void DriftReducedSystem::v_GenerateSummary(SU::SummaryList &s) {
  UnsteadySystem::v_GenerateSummary(s);

  std::stringstream tmpss;
  tmpss << "[" << m_d00 << "," << m_d11 << "," << m_d22 << "]";
  SU::AddSummaryItem(s, "Helmsolve coeffs.", tmpss.str());

  SU::AddSummaryItem(s, "Reference density", m_n_ref);
  SU::AddSummaryItem(s, "Density floor", m_n_floor_fac);

  SU::AddSummaryItem(s, "Riemann solver", m_riemann_solver_type);
  // Particle stuff
  SU::AddSummaryItem(s, "Num. part. substeps", m_num_part_substeps);
  SU::AddSummaryItem(s, "Part. output freq", m_num_write_particle_steps);
  tmpss = std::stringstream();
  tmpss << "[" << m_B[0] << "," << m_B[1] << "," << m_B[2] << "]";
  SU::AddSummaryItem(s, "B", tmpss.str());
  SU::AddSummaryItem(s, "|B|", m_Bmag);
}

/**
 * @brief Post-construction class initialisation.
 *
 * @param create_field if true, create a new field object and add it to
 * m_fields. Optional, defaults to true.
 */
void DriftReducedSystem::v_InitObject(bool create_field) {
  // If particle-coupling is enabled,
  if (this->m_particle_sys->m_num_particles > 0) {
    m_required_flds.push_back("ne_src");
  }

  AdvectionSystem::v_InitObject(create_field);

  // Ensure that the session file defines all required variables and that they
  // have the same order
  validate_fields();

  // Load parameters
  load_params();

  // Tell UnsteadySystem to only integrate a subset of fields in time
  // (Ignore fields that don't have a time derivative)
  m_intVariables.resize(m_int_fld_names.size());
  for (auto ii = 0; ii < m_int_fld_names.size(); ii++) {
    int var_idx = m_field_to_index.get_idx(m_int_fld_names[ii]);
    ASSERTL0(var_idx >= 0, "Setting time integration vars - GetIntFieldNames() "
                           "returned an invalid field name.");
    m_intVariables[ii] = var_idx;
  }

  // Since we are starting from a setup where each field is defined to be a
  // discontinuous field (and thus support DG), the first thing we do is to
  // recreate the phi field so that it is continuous, in order to support the
  // Poisson solve. Note that you can still perform a Poisson solve using a
  // discontinuous field, which is done via the hybridisable discontinuous
  // Galerkin (HDG) approach.
  int phi_idx = m_field_to_index.get_idx("phi");
  m_fields[phi_idx] = MemoryManager<MR::ContField>::AllocateSharedPtr(
      m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

  // Create storage for advection velocities, parallel velocity difference,ExB
  // drift velocity, E field
  int npts = GetNpoints();
  for (int i = 0; i < m_graph->GetSpaceDimension(); ++i) {
    m_adv_vel_elec[i] = Array<OneD, NekDouble>(npts);
    m_ExB_vel[i] = Array<OneD, NekDouble>(npts);
    m_E[i] = Array<OneD, NekDouble>(npts);
  }
  // Create storage for electron parallel velocities
  m_par_vel_elec = Array<OneD, NekDouble>(npts);

  // Type of advection class to be used. By default, we only support the
  // discontinuous projection, since this is the only approach we're
  // considering for this solver.
  ASSERTL0(m_projectionType == MR::eDiscontinuous,
           "Unsupported projection type: only discontinuous"
           " projection supported."); ////

  // Do not forwards transform initial condition.
  m_homoInitialFwd = false; ////

  // Define the normal velocity fields.
  // These are populated at each step (by reference) in calls to GetVnAdv()
  if (m_fields[0]->GetTrace()) {
    auto nTrace = GetTraceNpoints();
    m_norm_vel_elec = Array<OneD, NekDouble>(nTrace);
    m_norm_vel_vort = Array<OneD, NekDouble>(nTrace);
  }

  // Advection objects
  // Need one per advection velocity
  m_adv_elec = SU::GetAdvectionFactory().CreateInstance(m_adv_type, m_adv_type);
  m_adv_vort = SU::GetAdvectionFactory().CreateInstance(m_adv_type, m_adv_type);

  // Set callback functions to compute flux vectors
  m_adv_elec->SetFluxVector(&DriftReducedSystem::get_flux_vector_elec, this);
  m_adv_vort->SetFluxVector(&DriftReducedSystem::get_flux_vector_vort, this);

  // Create Riemann solvers (one per advection object) and set normal velocity
  // callback functions
  m_riemann_elec = SU::GetRiemannSolverFactory().CreateInstance(
      m_riemann_solver_type, m_session);
  m_riemann_elec->SetScalar("Vn", &DriftReducedSystem::get_adv_vel_norm_elec,
                            this);
  m_riemann_vort = SU::GetRiemannSolverFactory().CreateInstance(
      m_riemann_solver_type, m_session);
  m_riemann_vort->SetScalar("Vn", &DriftReducedSystem::get_adv_vel_norm_vort,
                            this);

  // Tell advection objects about the Riemann solvers and finish init
  m_adv_elec->SetRiemannSolver(m_riemann_elec);
  m_adv_elec->InitObject(m_session, m_fields);
  m_adv_vort->SetRiemannSolver(m_riemann_vort);
  m_adv_vort->InitObject(m_session, m_fields);

  // Bind projection function for time integration object
  m_ode.DefineProjection(&DriftReducedSystem::do_ode_projection, this);

  ASSERTL0(m_explicitAdvection,
           "This solver only supports explicit-in-time advection.");

  // Store DisContFieldSharedPtr casts of fields in a map, indexed by name, for
  // use in particle project,evaluate operations
  int idx = 0;
  for (auto &field_name : m_session->GetVariables()) {
    m_discont_fields[field_name] =
        std::dynamic_pointer_cast<MR::DisContField>(m_fields[idx]);
    idx++;
  }

  if (m_particle_sys->m_num_particles > 0) {
    // Set up object to project onto density source field
    int low_order_project;
    m_session->LoadParameter("low_order_project", low_order_project, 0);
    if (low_order_project) {
      ASSERTL0(
          m_discont_fields.count("ne_src_interp"),
          "Intermediate, lower order interpolation field not found in config.");
      m_particle_sys->setup_project(m_discont_fields["ne_src_interp"],
                                    m_discont_fields["ne_src"]);
    } else {
      m_particle_sys->setup_project(m_discont_fields["ne_src"]);
    }
  }

  // Set up object to evaluate density field
  m_particle_sys->setup_evaluate_ne(m_discont_fields["ne"]);
}

/**
 * @brief Override v_PostIntegrate to do particle output
 *
 * @param step Time step number
 */
bool DriftReducedSystem::v_PostIntegrate(int step) {
  // Writes a step of the particle trajectory.
  if (m_num_write_particle_steps > 0 &&
      (step % m_num_write_particle_steps) == 0) {
    m_particle_sys->write(step);
    m_particle_sys->write_source_fields();
  }
  return AdvectionSystem::v_PostIntegrate(step);
}

/**
 * @brief Override v_PreIntegrate to do particle system integration, projection
 * onto source terms.
 *
 * @param step Time step number
 */
bool DriftReducedSystem::v_PreIntegrate(int step) {
  if (m_particle_sys->m_num_particles > 0) {
    // Integrate the particle system to the requested time.
    m_particle_sys->integrate(m_time + m_timestep, m_part_timestep);
    // Project onto the source fields
    m_particle_sys->project_source_terms();
  }

  return AdvectionSystem::v_PreIntegrate(step);
}

/**
 * @brief Convenience function to zero a Nektar Array of 1D Arrays.
 *
 * @param out_arr Array of 1D arrays to be zeroed
 *
 */
void DriftReducedSystem::zero_out_array(
    Array<OneD, Array<OneD, NekDouble>> &out_arr) {
  for (auto ifld = 0; ifld < out_arr.size(); ifld++) {
    Vmath::Zero(out_arr[ifld].size(), out_arr[ifld], 1);
  }
}
} // namespace NESO::Solvers::H3LAPD

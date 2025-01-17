#ifndef __EQNSYS_BASE_H_
#define __EQNSYS_BASE_H_

#include <SolverUtils/EquationSystem.h>

#include <type_traits>

#include "nektar_interface/solver_base/partsys_base.hpp"
#include "nektar_interface/utilities.hpp"

namespace LU = Nektar::LibUtilities;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;

namespace NESO::Solvers {

// Fwd declare NEKEQNSYS, PARTSYS as classes
class NEKEQNSYS;
class PARTSYS;

/**
 * @brief Base class for Nektar++ equation systems, coupled to a NESO-Particles
 * particle system derived from NESO::Solvers::PartSysBase.
 */
template <typename NEKEQNSYS, typename PARTSYS>
class EqnSysBase : public NEKEQNSYS {
  // Equation system template param must derive from Nektar's EquationSystem
  // base class
  static_assert(std::is_base_of<SU::EquationSystem, NEKEQNSYS>(),
                "Template arg to EqnSysBase must derive from "
                "Nektar::SolverUtils::EquationSystem");
  // Particle system template param must derive from PartSysBase
  static_assert(
      std::is_base_of<PartSysBase, PARTSYS>(),
      "PARTSYS template arg to EqnSysBase must derive from PartSysBase");

protected:
  EqnSysBase(const LU::SessionReaderSharedPtr &session,
             const SD::MeshGraphSharedPtr &graph)
      : NEKEQNSYS(session, graph), field_to_index(session->GetVariables()),
        required_fld_names() {

    // If number of particles / number per cell was set in config; construct the
    // particle system
    int num_parts_per_cell, num_parts_tot;
    session->LoadParameter(PartSysBase::NUM_PARTS_TOT_STR, num_parts_tot, -1);
    session->LoadParameter(PartSysBase::NUM_PARTS_PER_CELL_STR,
                           num_parts_per_cell, -1);
    this->particles_enabled = num_parts_tot > 0 || num_parts_per_cell > 0;
    if (this->particles_enabled) {
      this->particle_sys = std::make_shared<PARTSYS>(session, graph);
    }
  }

  /// Field name => index mapper
  NESO::NektarFieldIndexMap field_to_index;

  /// Store mesh dims and number of quad points as member vars for convenience
  int n_dims;
  int n_pts;

  /// Particle system
  std::shared_ptr<PARTSYS> particle_sys;

  /// Flag identifying whether particles were enabled in the config file
  bool particles_enabled;

  /// List of field names required by the solver
  std::vector<std::string> required_fld_names;

  /// Placeholder for subclasses to override; called in v_InitObject()
  virtual void load_params(){};

  /// Hook to allow subclasses to run post-solve tasks at the end of v_DoSolve()
  virtual void post_solve(){};

  /**
   * @brief Assert that a named variable/field is at a particular index in the
   * sessions variable list.
   */
  void check_var_idx(const int &idx, const std::string var_name) {
    std::stringstream err;
    err << "Expected variable index " << idx << " to correspond to '"
        << var_name << "'. Check your session file.";
    NESOASSERT(this->m_session->GetVariable(idx).compare(var_name) == 0,
               err.str());
  }

  /**
   * @brief Check that all required fields are defined. All fields must have the
   * same number of quad points for now.
   */
  void validate_fields() {
    for (auto &fld_name : this->required_fld_names) {
      // Check field exists
      int idx = this->field_to_index.get_idx(fld_name);
      std::string err_msg = "Required field [" + fld_name + "] is not defined.";
      NESOASSERT(idx >= 0, err_msg.c_str());

      // Check fields all have the same number of quad points
      int n_pts = this->m_fields[idx]->GetNpoints();
      err_msg = "Expecting " + std::to_string(this->n_pts) +
                " quad points, but field '" + fld_name + "' has " +
                std::to_string(n_pts) +
                ". Check NUMMODES is the same for all required fields.";
      NESOASSERT(n_pts == this->n_pts, err_msg.c_str());
    }
  }

  /** @brief Write particle params to stdout on task 0. Ensures they appear just
   * after fluid params are written by nektar.
   *
   * @param dump_initial_conditions Write initial conditions to file?
   * default=true
   *
   * */
  virtual void v_DoInitialise(bool dump_initial_conditions) override {
    if (this->m_session->GetComm()->TreatAsRankZero() &&
        this->particles_enabled) {
      particle_sys->add_params_report();
    }
    NEKEQNSYS::v_DoInitialise(dump_initial_conditions);
  }

  /**
   * @brief Free particle system memory after solver loop has finished.
   * Prevent further overrides to guarantee that subclasses do the same.
   * Subclasses can override post_solve() to add additional tasks.
   */
  virtual void v_DoSolve() override final {
    NEKEQNSYS::v_DoSolve();
    if (this->particle_sys) {
      this->particle_sys->free();
    }
    post_solve();
  }

  /**
   * @brief Initialise the equation system, then check required fields are set
   * and load parameters.
   */
  virtual void v_InitObject(bool create_fields) override {
    NEKEQNSYS::v_InitObject(create_fields);

    this->n_dims = NEKEQNSYS::m_graph->GetMeshDimension();
    this->n_pts = NEKEQNSYS::m_fields[0]->GetNpoints();

    // Ensure that the session file defines all required variables
    validate_fields();

    // Load parameters
    load_params();
  }

  /**
   * @brief Utility function to zero a Nektar++ array of arrays
   *
   * @param arr Array-of-arrays to zero
   */
  void zero_array_of_arrays(Array<OneD, Array<OneD, NekDouble>> &arr) {
    for (auto ii = 0; ii < arr.size(); ii++) {
      Vmath::Zero(arr[ii].size(), arr[ii], 1);
    }
  }
};

} // namespace NESO::Solvers
#endif
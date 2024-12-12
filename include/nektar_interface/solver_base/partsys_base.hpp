#ifndef __PARTSYS_BASE_H_
#define __PARTSYS_BASE_H_

#include <SolverUtils/EquationSystem.h>
#include <mpi.h>
#include <nektar_interface/geometry_transport/halo_extension.hpp>
#include <nektar_interface/particle_interface.hpp>
#include <neso_particles.hpp>
#include <type_traits>

namespace LU = Nektar::LibUtilities;
namespace SD = Nektar::SpatialDomains;

namespace NESO::Particles {

/// Struct used to set common options for particle systems
struct PartSysOptions {
  int extend_halos_offset = 0;
};

class PartSysBase;

/// Nektar style factory method
typedef ParticleSystemSharedPtr std::shared_ptr<PartSysBase>;
typedef ParticleSystemFactory
    LU::NekFactory<std::string, PartSysBase, const LU::SessionReaderSharedPtr,
                   const SD::MeshGraphSharedPtr>;
ParticleSystemFactory &GetParticleSystemFactory();

class PartSysBase {

public:
  static std::string className;

  // Some parameter names used in solver config files
  inline static const std::string NUM_PARTS_TOT_STR = "num_particles_total";
  inline static const std::string NUM_PARTS_PER_CELL_STR =
      "num_particles_per_cell";
  inline static const std::string PART_OUTPUT_FREQ_STR = "particle_output_freq";

  /// Total number of particles in simulation
  int64_t num_parts_tot;

  /// NESO-Particles ParticleSpec;
  ParticleSpec particle_spec;
  /// NESO-Particles ParticleGroup
  ParticleGroupSharedPtr particle_group;

  /// Compute target
  SYCLTargetSharedPtr sycl_target;

  /// @brief Report particle parameter values (writes to stdout).
  void add_params_report();

  /// @brief Clear up memory related to the particle system
  void free();

  /**
   * @brief Check whether particle output is scheduled for \p step.
   *
   * @param step
   * @returns true if \p step is a scheduled output step, according to the
   * frequency read from the config file, false otherwise
   */
  bool is_output_step(int step);

  /**
   *  @brief Write particle properties to an output file.
   *  @param step Time step number.
   */
  void write(const int step);

  /// @brief Reads the particle tag from xml document
  virtual void ReadParticles();

  /// @brief  Reads parameters related to particles
  /// @param particles
  void ReadParameters(TiXmlElement *particles);

  /// @brief Reads the list of species defined under particles
  /// @param particles
  void ReadSpecies(TiXmlElement *particles);

  /// @brief Reads the particle boundary conditions
  /// @params particles
  void ReadBoundary(TiXmlElement *particles);

  // Make this pure virtual?
  /// @brief Instantiates the particle spec
  virtual void InitSpec();

  /// @brief Instantiates the particle system object, including the
  /// particle_group.  Delayed until after spec is determined from reading xml
  virtual void InitObject();

protected:
  /**
   * @brief Protected constructor to prohibit direct instantiation.
   *  @param session Nektar++ session to use for parameters and simulation
   * specification.
   *  @param graph Nektar++ MeshGraph on which particles exist.
   *  @param comm (optional) MPI communicator to use - default MPI_COMM_WORLD.
   *
   */
  PartSysBase(const LU::SessionReaderSharedPtr session,
              const SD::MeshGraphSharedPtr graph,
              MPI_Comm comm = MPI_COMM_WORLD,
              PartSysOptions options = PartSysOptions());

  /// Object used to map to/from nektar geometry ids to 0,N-1
  std::shared_ptr<CellIDTranslation> cell_id_translation;
  /// MPI communicator
  MPI_Comm comm;
  /// NESO-Particles domain.
  DomainSharedPtr domain;
  /// Pointer to Nektar Meshgraph object
  SD::MeshGraphSharedPtr graph;
  /// HDF5 output file
  std::shared_ptr<H5Part> h5part;
  /// Number of spatial dimensions being used
  const int ndim;
  /// Mapping instance to map particles into nektar++ elements.
  std::shared_ptr<NektarGraphLocalMapper> nektar_graph_local_mapper;
  /// Options struct
  PartSysOptions options;
  /// HMesh instance that allows particles to move over nektar++ meshes.
  ParticleMeshInterfaceSharedPtr particle_mesh_interface;
  /// Pointer to Session object
  LU::SessionReaderSharedPtr session;

  ParameterMap m_parameters;
  InterpreterSharedPtr m_interpreter;

  /**
   * @brief Set up per-step particle output
   *  @param fname Output filename. Default is 'particle_trajectory.h5part'.
   *  @param args Remaining arguments (variable length) should be sym instances
   *  indicating which ParticleDats are to be written.
   */
  template <typename... T> void init_output(std::string fname, T... args) {
    if (this->h5part) {
      if (this->sycl_target->comm_pair.rank_parent == 0) {
        nprint("Ignoring (duplicate?) call to init_output().");
      }
    } else {
      // Create H5Part instance
      this->h5part =
          std::make_shared<H5Part>(fname, this->particle_group, args...);
    }
  }

  /// @brief Read some parameters associated with all particle systems.
  void read_params();

  /**
   * @brief Store particle param values in a map.Values are reported later via
   * add_params_report()
   * @param label Label to attach to the parameter
   * @param value Value of the parameter that was set
   */
  template <typename T> void report_param(std::string label, T val) {
    // form stringstream and store string value in private map
    std::stringstream ss;
    ss << val;
    param_vals_to_report[label] = ss.str();
  }

private:
  /// Output frequency read from config file
  int output_freq;
  /**
   * Map containing parameter name,value pairs to be written to stdout when
   * the nektar equation system is initialised. Populated with report_param().
   */
  std::map<std::string, std::string> param_vals_to_report;
};

} // namespace NESO::Particles
#endif
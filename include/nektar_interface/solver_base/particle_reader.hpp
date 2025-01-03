#ifndef __PARTICLE_READER_H_
#define __PARTICLE_READER_H_

#include <LibUtilities/BasicUtils/SessionReader.h>

class TiXmlElement;
class TiXmlDocument;
class TiXmlHandle;

using namespace Nektar;
namespace LU = Nektar::LibUtilities;

namespace NESO::Particles {

typedef std::map<std::string, std::string> SpeciesParamMap;
typedef std::pair<SpeciesParamMap, LU::FunctionMap> SpeciesMap;
typedef std::map<std::string, SpeciesMap> SpeciesMapList;

enum class ParticleBoundaryConditionType {
  ePeriodic,
  eReflective,
  eNotDefined
};

typedef std::map<std::string, ParticleBoundaryConditionType>
    SpeciesBoundaryList;
typedef std::map<int, SpeciesBoundaryList> ParticleBoundaryList;

class ParticleReader;
typedef std::shared_ptr<ParticleReader> ParticleReaderSharedPtr;

class ParticleReader {
public:
  ParticleReader(const LU::SessionReaderSharedPtr session)
      : m_session(session) {};

  /// @brief Reads the particle tag from xml document
  void ReadParticles();
  /// @brief Reads info related to particles
  void ReadInfo();
  const std::string &GetInfo(const std::string &pName) const;

  /// @brief  Reads parameters related to particles
  /// @param particles
  void ReadParameters(TiXmlElement *particles);

  /// @brief  Reads functions related to a species (e.g. Initial Conditions)
  /// @param particles
  /// @param functions
  void ReadSpeciesFunctions(TiXmlElement *particles,
                            LU::FunctionMap &functions);

  /// @brief Reads the list of species defined under particles
  /// @param particles
  void ReadSpecies(TiXmlElement *particles);

  /// @brief Reads the particle boundary conditions
  /// @params particles
  void ReadBoundary(TiXmlElement *particles);

  /// Checks if a parameter is specified in the XML document.
  bool DefinesParameter(const std::string &name) const;
  /// Returns the value of the specified parameter.
  const NekDouble &GetParameter(const std::string &pName) const;

  /// Load an integer parameter
  void LoadParameter(const std::string &name, int &var) const;
  /// Load an size_t parameter
  void LoadParameter(const std::string &name, size_t &var) const;
  /// Check for and load an integer parameter.
  void LoadParameter(const std::string &name, int &var, const int &def) const;
  /// Check for and load an size_t parameter.
  void LoadParameter(const std::string &name, size_t &var,
                     const size_t &def) const;
  /// Load a double precision parameter
  void LoadParameter(const std::string &name, NekDouble &var) const;
  /// Check for and load a double-precision parameter.
  void LoadParameter(const std::string &name, NekDouble &var,
                     const NekDouble &def) const;

private:
  LU::SessionReaderSharedPtr m_session;
  /// Map of particle info (e.g. Particle System name)
  std::map<std::string, std::string> m_particleInfo;
  // Map of specied
  SpeciesMapList m_species;
  LU::ParameterMap m_parameters;
  LU::InterpreterSharedPtr m_interpreter;
  /// Functions.
  LU::FunctionMap m_functions;

  ParticleBoundaryList m_boundaryConditions;

  void ParseEquals(const std::string &line, std::string &lhs, std::string &rhs);
};

} // namespace NESO::Particles
#endif
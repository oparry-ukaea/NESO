///////////////////////////////////////////////////////////////////////////////
//
// File ReducedH3LAPDSystem.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Header for the Hermes-3 LAPD equation system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ReducedH3LAPDSystem_H
#define ReducedH3LAPDSystem_H

#include "nektar_interface/utilities.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar {

class ReducedH3LAPDSystem : virtual public SolverUtils::AdvectionSystem {
public:
  friend class MemoryManager<ReducedH3LAPDSystem>;

  /// Name of class.
  static std::string className;

  /// Creates an instance of this class.
  static SolverUtils::EquationSystemSharedPtr
  create(const LibUtilities::SessionReaderSharedPtr &pSession,
         const SpatialDomains::MeshGraphSharedPtr &pGraph) {
    SolverUtils::EquationSystemSharedPtr p =
        MemoryManager<ReducedH3LAPDSystem>::AllocateSharedPtr(pSession, pGraph);
    p->InitObject();
    return p;
  }

  /// Default destructor.
  virtual ~ReducedH3LAPDSystem() = default;

protected:
  /// Protected constructor. Since we use a factory pattern, objects should be
  /// constructed via the SolverUtils::EquationSystem factory.
  ReducedH3LAPDSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::MeshGraphSharedPtr &pGraph);

  // Field name => index mapper
  NESO::NektarFieldIndexMap m_field_to_index;
  // Forcing/source terms
  std::vector<SolverUtils::ForcingSharedPtr> m_forcing;
  // List of field names required by the solver
  std::vector<std::string> m_required_flds;

  void AddAdvTerms(std::vector<std::string> field_names,
                   const SolverUtils::AdvectionSharedPtr advObj,
                   const Array<OneD, Array<OneD, NekDouble>> &vAdv,
                   const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                   Array<OneD, Array<OneD, NekDouble>> &outarray,
                   const NekDouble time);

  void AddCollisionAndPolDriftTerms(
      const Array<OneD, const Array<OneD, NekDouble>> &inarray,
      Array<OneD, Array<OneD, NekDouble>> &outarray);
  void AddEParTerms(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                    Array<OneD, Array<OneD, NekDouble>> &outarray);
  void AddGradPTerms(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);
  void AddDivvParTerm(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                      Array<OneD, Array<OneD, NekDouble>> &outarray);
  void CalcCollisionFreqs(const Array<OneD, NekDouble> &ne,
                          Array<OneD, NekDouble> &coeffs);
  void CalcCoulombLogarithm(const Array<OneD, NekDouble> &ne,
                            Array<OneD, NekDouble> &LogLambda);
  void
  CalcEAndAdvVels(const Array<OneD, const Array<OneD, NekDouble>> &inarray);
  void CalcAdvNormalVels();
  void DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time);

  void ExplicitTimeInt(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time);
  void GetFluxVectorDiff(
      const Array<OneD, Array<OneD, NekDouble>> &inarray,
      const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
      Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor);

  void
  GetFluxVectorElec(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);
  void
  GetFluxVectorVort(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

  Array<OneD, NekDouble> &
  GetVnAdv(Array<OneD, NekDouble> &traceVn,
           const Array<OneD, Array<OneD, NekDouble>> &vAdv);

  Array<OneD, NekDouble> &GetVnAdvElec();
  Array<OneD, NekDouble> &GetVnAdvVort();

  void LoadParams();

  void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &inarray);

  void ValidateFieldList();

  virtual void v_InitObject(bool DeclareField) override;

private:
  NekDouble m_d22;
  // Advection type
  std::string m_advType;
  // Magnetic field vector
  std::vector<NekDouble> m_B;
  // Magnitude of the magnetic field
  NekDouble m_Bmag;
  // Normalised magnetic field vector
  std::vector<NekDouble> m_b_unit;
  // Charge unit
  NekDouble m_charge_e;
  // Ion mass;
  NekDouble m_md;
  // Electron mass;
  NekDouble m_me;
  // Reference number density
  NekDouble m_nRef;
  // Riemann solver type (used for all advection terms)
  std::string m_RiemSolvType;
  // Ion temperature in eV
  NekDouble m_Td;
  // Electron temperature in eV
  NekDouble m_Te;
  //---------------------------------------------------------------------------
  // Factors used in collision coeff calculation
  // Density-independent part of the Coulomb logarithm; read from config
  NekDouble m_coulomb_log_const;
  // Pre-factor used when calculating collision frequencies; read from config
  NekDouble m_nu_ei_const;
  // Factor to convert densities (back) to SI; used in Coulomb logarithm calc
  NekDouble m_n_to_SI;
  //---------------------------------------------------------------------------
  // Advection objects
  SolverUtils::AdvectionSharedPtr m_advElec;
  SolverUtils::AdvectionSharedPtr m_advVort;
  // Storage for Electric field
  Array<OneD, Array<OneD, NekDouble>> m_E;
  // Riemann solver objects
  SolverUtils::RiemannSolverSharedPtr m_riemannSolverElec;
  SolverUtils::RiemannSolverSharedPtr m_riemannSolverVort;
  // Storage for advection velocities dotted with element_edge_normals
  Array<OneD, NekDouble> m_traceVnElec;
  Array<OneD, NekDouble> m_traceVnVort;
  // Storage for electron advection velocity
  Array<OneD, Array<OneD, NekDouble>> m_vAdvElec;
  // Storage for ExB drift velocity
  Array<OneD, Array<OneD, NekDouble>> m_vExB;
  // Storage for electron perpendicular velocity
  Array<OneD, NekDouble> m_vPerpElec;
  //---------------------------------------------------------------------------
  // Debugging
  void PrintArrVals(const Array<OneD, NekDouble> &arr, int num, int stride = 1,
                    std::string label = "", bool all_tasks = false);
  void PrintArrSize(Array<OneD, NekDouble> &arr, std::string label = "",
                    bool all_tasks = false);
};

} // namespace Nektar
#endif
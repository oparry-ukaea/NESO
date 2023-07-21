///////////////////////////////////////////////////////////////////////////////
//
// File ReducedH3LAPDSystem.cpp
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
// Description: Hermes-3 LAPD equation system
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "ReducedH3LAPDSystem.h"

namespace Nektar {
std::string ReducedH3LAPDSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "reducedLAPD", ReducedH3LAPDSystem::create,
        "Reduced Hermes-3 LAPD equation system");

ReducedH3LAPDSystem::ReducedH3LAPDSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      m_field_to_index(pSession->GetVariables()),
      m_vAdvElec(pGraph->GetSpaceDimension()),
      m_vExB(pGraph->GetSpaceDimension()), m_E(pGraph->GetSpaceDimension()) {
  m_required_flds = {"ne", "w", "phi"};
}

/**
 * @brief Perform projection into correct polynomial space.
 *
 * This routine projects the @p inarray input and ensures the @p outarray
 * output lives in the correct space. Since we are hard-coding DG, this
 * corresponds to a simple copy from in to out, since no elemental
 * connectivity is required and the output of the RHS function is
 * polynomial.
 */
void ReducedH3LAPDSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {
  int nvariables = inarray.size();
  int npoints = inarray[0].size();
  SetBoundaryConditions(time);

  for (int i = 0; i < nvariables; ++i) {
    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
  }
}

void ReducedH3LAPDSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {
  int nPts = GetNpoints();

  // Set advection vel
  Vmath::Fill(nPts, 0.0, m_vAdvElec[0], 1);
  Vmath::Fill(nPts, 0.0, m_vAdvElec[1], 1);
  Vmath::Fill(nPts, 0.05, m_vAdvElec[2], 1);

  Array<OneD, Array<OneD, NekDouble>> tmp_out(nPts);
  tmp_out[0] = Array<OneD, NekDouble>(nPts, 0.0);
  m_advElec->Advect(1, m_fields, m_vAdvElec, inarray, tmp_out, time);
  Vmath::Vsub(nPts, outarray[0], 1, tmp_out[0], 1, outarray[0], 1);
}

/**
 * @brief Compute the flux vector for advecting the electron density
 *
 * @param physfield   Array of Fields ptrs
 * @param flux        Resulting flux array
 */
void ReducedH3LAPDSystem::GetFluxVectorElec(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {

  ASSERTL1(flux[0].size() == vAdv.size(),
           "Dimension of flux array and advection velocity array do not match");
  int nq = physfield[0].size();

  for (auto ifld = 0; ifld < flux.size(); ++ifld) {
    for (auto idim = 0; idim < flux[0].size(); ++idim) {
      Vmath::Vmul(nq, physfield[ifld], 1, m_vAdvElec[idim], 1, flux[ifld][idim],
                  1);
    }
  }
}

/**
 * @brief Compute the normal advection velocity for the electron density
 */
Array<OneD, NekDouble> &ReducedH3LAPDSystem::GetVnAdvElec() {
  // Number of trace (interface) points
  int nTracePts = GetTraceNpoints();
  // Auxiliary variable to compute normal velocities
  Array<OneD, NekDouble> tmp(nTracePts);

  // Zero previous values
  Vmath::Zero(nTracePts, m_traceVnElec, 1);

  //  Compute dot product of advection velocity with the trace normals and store
  for (int i = 0; i < m_vAdvElec.size(); ++i) {
    m_fields[0]->ExtractTracePhys(m_vAdvElec[i], tmp);
    Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, m_traceVnElec, 1,
                 m_traceVnElec, 1);
  }
  return m_traceVnElec;
}

void ReducedH3LAPDSystem::LoadParams() {
  // Type of advection to use -- in theory we also support flux reconstruction
  // for quad-based meshes, or you can use a standard convective term if you
  // were fully continuous in space. Default is DG.
  m_session->LoadSolverInfo("AdvectionType", m_advType, "WeakDG");

  // ***Assumes field aligned with z-axis***
  // Magnetic field strength. Fix B = [0, 0, Bxy] for now
  m_B = std::vector<NekDouble>(m_graph->GetSpaceDimension(), 0);
  m_session->LoadParameter("Bxy", m_B[2], 0.1);

  // Charge
  m_session->LoadParameter("e", m_charge_e, 1.0);

  // Ion mass
  m_session->LoadParameter("md", m_md, 2.0);

  // Electron mass - default val is multiplied by 60 to improve convergence
  m_session->LoadParameter("me", m_me, 60. / 1836);

  // Reference number density
  m_session->LoadParameter("nRef", m_nRef, 1.0);

  // Electron temperature in eV
  m_session->LoadParameter("Te", m_Te, 5.0);

  // Ion temperature in eV
  m_session->LoadParameter("Td", m_Td, 0.1);

  // Type of Riemann solver to use. Default = "Upwind"
  m_session->LoadSolverInfo("UpwindType", m_RiemSolvType, "Upwind");

  // d22
  m_session->LoadParameter("d22", m_d22, 0.0);
}

void ReducedH3LAPDSystem::PrintArrSize(Array<OneD, NekDouble> &arr,
                                       std::string label, bool all_tasks) {
  if (m_session->GetComm()->TreatAsRankZero() || all_tasks) {
    if (!label.empty()) {
      std::cout << label << " ";
    }
    std::cout << "size = " << arr.size() << std::endl;
  }
}

void ReducedH3LAPDSystem::PrintArrVals(const Array<OneD, NekDouble> &arr,
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
 * Check all required fields are defined
 */
void ReducedH3LAPDSystem::ValidateFieldList() {
  for (auto &fld_name : m_required_flds) {
    ASSERTL0(m_field_to_index.get_idx(fld_name) >= 0,
             "Required field [" + fld_name + "] is not defined.");
  }
}

/**
 * @brief Initialization object for ReducedH3LAPDSystem class.
 */
void ReducedH3LAPDSystem::v_InitObject(bool DeclareField) {
  //  Ensure that the session file defines all required variables
  ValidateFieldList();

  // Load parameters
  LoadParams();

  // Compute some properties derived from params
  m_Bmag = std::sqrt(m_B[0] * m_B[0] + m_B[1] * m_B[1] + m_B[2] * m_B[2]);
  m_b_unit = std::vector<NekDouble>(m_graph->GetSpaceDimension());
  for (auto idim = 0; idim < m_b_unit.size(); idim++) {
    m_b_unit[idim] = m_B[idim] / m_Bmag;
  }

  AdvectionSystem::v_InitObject(DeclareField);

  // Tell UnsteadySystem to only integrate a subset of fields in time
  // (i.e. electron density, parallel momenta and vorticity), ignoring
  // electrostatic potential, since this doesn't have a time-derivative.
  int nVar = m_fields.size();
  m_intVariables.resize(nVar - 1);
  for (int i = 0; i < nVar - 1; ++i) {
    m_intVariables[i] = i;
  }

  // Since we are starting from a setup where each field is defined to be a
  // discontinuous field (and thus support DG), the first thing we do is to
  // recreate the phi field so that it is continuous, in order to support the
  // Poisson solve. Note that you can still perform a Poisson solve using a
  // discontinuous field, which is done via the hybridisable discontinuous
  // Galerkin (HDG) approach.
  int phi_idx = m_field_to_index.get_idx("phi");
  m_fields[phi_idx] = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
      m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

  // Create storage for total advection velocity, ExB drift velocity, E field
  int nPts = GetNpoints();
  for (int i = 0; i < m_graph->GetSpaceDimension(); ++i) {
    m_vAdvElec[i] = Array<OneD, NekDouble>(nPts);
    m_vExB[i] = Array<OneD, NekDouble>(nPts);
    m_E[i] = Array<OneD, NekDouble>(nPts);
  }

  // Type of advection class to be used. By default, we only support the
  // discontinuous projection, since this is the only approach we're
  // considering for this solver.
  ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
           "Unsupported projection type: only discontinuous"
           " projection supported."); ////

  // Do not forwards transform initial condition.
  m_homoInitialFwd = false; ////

  // Define the normal velocity fields.
  if (m_fields[0]->GetTrace()) {
    auto nTrace = GetTraceNpoints();
    m_traceVnElec = Array<OneD, NekDouble>(nTrace);
  }

  // Advection objects
  // Need one per advection velocity - for (ne,Ge), Gd and w
  m_advElec =
      SolverUtils::GetAdvectionFactory().CreateInstance(m_advType, m_advType);

  // Set callback functions to compute flux vectors
  m_advElec->SetFluxVector(&ReducedH3LAPDSystem::GetFluxVectorElec, this);

  // Create Riemann solvers (one per advection object) and set normal velocity
  // callback functions
  m_riemannSolverElec = SolverUtils::GetRiemannSolverFactory().CreateInstance(
      m_RiemSolvType, m_session);
  m_riemannSolverElec->SetScalar("Vn", &ReducedH3LAPDSystem::GetVnAdvElec,
                                 this);

  // Tell advection objects about the Riemann solvers and finish init
  m_advElec->SetRiemannSolver(m_riemannSolverElec);
  m_advElec->InitObject(m_session, m_fields);

  // The m_ode object defines the timestepping to be used, and lives in
  // the SolverUtils::UnsteadySystem class. For explicit solvers, you need
  // to supply a right-hand side function, and a projection function
  // (e.g. for continuous Galerkin this would be an assembly-type
  // operation to ensure C^0 connectivity). These are done again through
  // callbacks.
  m_ode.DefineOdeRhs(&ReducedH3LAPDSystem::ExplicitTimeInt, this);
  m_ode.DefineProjection(&ReducedH3LAPDSystem::DoOdeProjection, this);

  ASSERTL0(m_explicitAdvection,
           "This solver only supports explicit-in-time advection.");
}

} // namespace Nektar

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
// Description: Reduced Hermes-3 LAPD equation system
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "CheckPhiH3LAPDSystem.h"
#include <filesystem>

namespace Nektar {
std::string CheckPhiH3LAPDSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "checkPhiLAPD", CheckPhiH3LAPDSystem::create,
        "checkPhi Hermes-3 LAPD equation system");

CheckPhiH3LAPDSystem::CheckPhiH3LAPDSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      H3LAPDSystem(pSession, pGraph) {
  m_required_flds = {"ne", "w", "phi"};
}

void CheckPhiH3LAPDSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {

  // Zero outarray
  for (auto ifld = 0; ifld < outarray.size(); ifld++) {
    Vmath::Zero(outarray[ifld].size(), outarray[ifld], 1);
  }

  // Solver for electrostatic potential.
  SolvePhi(inarray);

  // Write phi to a separate vtu
  std::filesystem::path phi_output_fname = std::filesystem::path("phi.vtu");
  int phi_idx = m_field_to_index.get_idx("phi");
  if (!std::filesystem::exists(phi_output_fname)) {
    NESO::write_vtu(m_fields[phi_idx], phi_output_fname.string(), "phi");
  }

  // Write phi diff to a separate vtu
  std::filesystem::path phi_diff_output_fname =
      std::filesystem::path("phi_diff.vtu");
  int nPts = GetNpoints();
  int phi_diff_idx = m_field_to_index.get_idx("phi_diff");
  int phi_sln_idx = m_field_to_index.get_idx("phi_sln");
  if (!std::filesystem::exists(phi_diff_output_fname)) {
    Vmath::Vsub(nPts, m_fields[phi_sln_idx]->GetPhys(), 1,
                m_fields[phi_idx]->GetPhys(), 1,
                m_fields[phi_diff_idx]->UpdatePhys(), 1);
    m_fields[phi_diff_idx]->FwdTrans(m_fields[phi_diff_idx]->GetPhys(),
                                     m_fields[phi_diff_idx]->UpdateCoeffs());
    NESO::write_vtu(m_fields[phi_diff_idx], phi_diff_output_fname.string(),
                    "phi_diff");
  }
}

void CheckPhiH3LAPDSystem::SolvePhi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray) {
  int nPts = GetNpoints();

  // Field indices
  int ne_idx = m_field_to_index.get_idx("ne");
  int phi_idx = m_field_to_index.get_idx("phi");
  int w_idx = m_field_to_index.get_idx("w");

  // Set up variable coefficients
  // ***Assumes field aligned with z-axis***
  StdRegions::VarCoeffMap varcoeffs;
  varcoeffs[StdRegions::eVarCoeffD00] = Array<OneD, NekDouble>(nPts, 1.0);
  varcoeffs[StdRegions::eVarCoeffD01] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD02] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD11] = Array<OneD, NekDouble>(nPts, 1.0);
  varcoeffs[StdRegions::eVarCoeffD12] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD22] = Array<OneD, NekDouble>(nPts, 0.0);

  // Set up factors for electrostatic potential solve. We support a generic
  // Helmholtz solve of the form (\nabla^2 - \lambda) u = f, so this sets
  // \lambda to zero.
  StdRegions::ConstFactorMap factors;
  factors[StdRegions::eFactorLambda] = 0.0;

  // Solve for phi. Output of this routine is in coefficient (spectral)
  // space, so backwards transform to physical space since we'll need that
  // for the advection step & computing drift velocity.
  m_fields[phi_idx]->HelmSolve(
      inarray[w_idx], m_fields[phi_idx]->UpdateCoeffs(), factors, varcoeffs);
  m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(),
                              m_fields[phi_idx]->UpdatePhys());
}

} // namespace Nektar

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
      H3LAPDSystem(pSession, pGraph) {
  m_required_flds = {"ne", "w", "phi"};
}

void ReducedH3LAPDSystem::AddDivvParTerm(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray) {
  int ne_idx = m_field_to_index.get_idx("ne");
  int w_idx = m_field_to_index.get_idx("w");
  int nPts = GetNpoints();
  Array<OneD, NekDouble> div_nvpar_term(nPts);
  Vmath::Vmul(nPts, inarray[ne_idx], 1, m_vParElec, 1, div_nvpar_term, 1);
  m_fields[ne_idx]->PhysDeriv(2, div_nvpar_term, div_nvpar_term);
  Vmath::Vsub(nPts, outarray[w_idx], 1, div_nvpar_term, 1, outarray[w_idx], 1);
}

/**
 * @brief Compute E = \f$ -\nabla\phi\f$, \f$ v_{E\times B}\f$ and the advection
 * velocities used in the ne/Ge, Gd equations.
 * @param inarray array of field physvals
 */
void ReducedH3LAPDSystem::CalcEAndAdvVels(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray) {
  int phi_idx = m_field_to_index.get_idx("phi");
  int nPts = GetNpoints();
  m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(), m_E[0], m_E[1],
                               m_E[2]);
  Vmath::Neg(nPts, m_E[0], 1);
  Vmath::Neg(nPts, m_E[1], 1);
  Vmath::Neg(nPts, m_E[2], 1);

  // v_ExB = Evec x Bvec / B^2
  Vmath::Svtsvtp(nPts, m_B[2] / m_Bmag / m_Bmag, m_E[1], 1,
                 -m_B[1] / m_Bmag / m_Bmag, m_E[2], 1, m_vExB[0], 1);
  Vmath::Svtsvtp(nPts, m_B[0] / m_Bmag / m_Bmag, m_E[2], 1,
                 -m_B[2] / m_Bmag / m_Bmag, m_E[0], 1, m_vExB[1], 1);
  Vmath::Svtsvtp(nPts, m_B[1] / m_Bmag / m_Bmag, m_E[0], 1,
                 -m_B[0] / m_Bmag / m_Bmag, m_E[1], 1, m_vExB[2], 1);

  // vpar
  int ne_idx = m_field_to_index.get_idx("ne");
  Array<OneD, NekDouble> tmpx(nPts), tmpy(nPts), tmpz(nPts);
  m_fields[ne_idx]->GetCoords(tmpx, tmpy, tmpz);
  m_session->GetFunction("vpar", ne_idx)
      ->Evaluate(tmpx, tmpy, tmpz, m_vParElec);
  for (auto iDim = 0; iDim < m_graph->GetSpaceDimension(); iDim++) {
    Vmath::Svtvp(nPts, m_b_unit[iDim], m_vParElec, 1, m_vExB[iDim], 1,
                 m_vAdvElec[iDim], 1);
  }
}

void ReducedH3LAPDSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {

  for (auto &var : {"ne", "Gd", "Ge", "w"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < inarray[fidx].size(); ii++) {
      if (!std::isfinite(inarray[fidx][ii])) {
        std::cout << "NaN in field " << var << ", aborting." << std::endl;
        exit(1);
      }
    }
  }

  // Zero outarrays
  for (auto ifld = 0; ifld < outarray.size(); ifld++) {
    Vmath::Zero(outarray[ifld].size(), outarray[ifld], 1);
  }

  // Solver for electrostatic potential.
  SolvePhi(inarray);

  // Calculate electric field from Phi, and compute v_ExB
  CalcEAndAdvVels(inarray);

  // Add ne advection term to outarray
  AddAdvTerms({"ne"}, m_advElec, m_vAdvElec, inarray, outarray, time);
  // Add w advection term to outarray
  AddAdvTerms({"w"}, m_advVort, m_vExB, inarray, outarray, time);

  AddDivvParTerm(inarray, outarray);
}

void ReducedH3LAPDSystem::LoadParams() {
  m_session->LoadSolverInfo("AdvectionType", m_advType, "WeakDG");

  // Magnetic field strength. Fix B = [0, 0, Bxy] for now
  m_B = std::vector<NekDouble>(m_graph->GetSpaceDimension(), 0);
  m_session->LoadParameter("Bxy", m_B[2], 0.1);

  // d22
  m_session->LoadParameter("d22", m_d22, 0.0);

  // Reference number density
  m_session->LoadParameter("nRef", m_nRef, 1.0);

  // Type of Riemann solver to use. Default = "Upwind"
  m_session->LoadSolverInfo("UpwindType", m_RiemSolvType, "Upwind");
}

void ReducedH3LAPDSystem::SolvePhi(
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
  varcoeffs[StdRegions::eVarCoeffD22] = Array<OneD, NekDouble>(nPts, m_d22);

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

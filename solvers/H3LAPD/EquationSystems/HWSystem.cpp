///////////////////////////////////////////////////////////////////////////////
//
// File HWSystem.cpp
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
// Description: Hasegawa-Waketani equation system as an intermediate step
// towards the full H3-LAPD problem.  Implemented by Ed Threlfall on 28/7/2023
// following email from Ben Dudson. System is basically advection + parallel
// Ohm's law only one parameter: \omega_{ce} / \nu_{ei} which is `large' ...
// evolves ne, w, phi only, no momenta, no ions
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "HWSystem.h"

namespace Nektar {
std::string HWSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "HWLAPD", HWSystem::create,
        "Hasegawa-Waketani equation system as an intermediate step towards the "
        "full H3-LAPD problem");

HWSystem::HWSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      H3LAPDSystem(pSession, pGraph) {
  m_required_flds = {"ne", "Gd", "Ge", "w", "phi"};
}

/**
 * Override CalcEAndAdvVels in order to set m_vAdvElec = m_vExB
 */
void HWSystem::CalcEAndAdvVels(
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

  // Set electron advection velocity = vExB
  // (m_vAdvIons isn't used, so isn't calculated here)
  for (auto iDim = 0; iDim < m_graph->GetSpaceDimension(); iDim++) {
    Vmath::Vcopy(nPts, m_vExB[iDim], 1, m_vAdvElec[iDim], 1);
  }
}

void HWSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {

  // Check inarray for NaNs
  for (auto &var : {"ne", "w"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < inarray[fidx].size(); ii++) {
      if (!std::isfinite(inarray[fidx][ii])) {
        std::cout << "NaN in field " << var << ", aborting." << std::endl;
        exit(1);
      }
    }
  }

  ZeroOutArray(outarray);

  // Solver for electrostatic potential.
  SolvePhi(inarray);

  // Calculate E from Phi, as well as corresponding advection velocities
  CalcEAndAdvVels(inarray);

  // Get field indices
  int nPts = GetNpoints();
  int ne_idx = m_field_to_index.get_idx("ne");
  int Ge_idx = m_field_to_index.get_idx("Ge");
  int Gd_idx = m_field_to_index.get_idx("Gd");
  int phi_idx = m_field_to_index.get_idx("phi");
  int w_idx = m_field_to_index.get_idx("w");

  // Advect both ne and w using ExB velocity
  AddAdvTerms({"ne"}, m_advElec, m_vExB, inarray, outarray, time);
  AddAdvTerms({"w"}, m_advVort, m_vExB, inarray, outarray, time);

  // Terms from 3D Hasegawa-Wakatani
  // Compute d^2 n/dz^2
  Array<OneD, NekDouble> first_deriv(nPts), HWterm_ne(nPts);
  m_fields[ne_idx]->PhysDeriv(2, m_fields[ne_idx]->GetPhys(), first_deriv);
  m_fields[ne_idx]->PhysDeriv(2, first_deriv, HWterm_ne);
  // Compute d^2 phi / dz^2
  Array<OneD, NekDouble> HWterm_phi(nPts);
  m_fields[ne_idx]->PhysDeriv(2, m_fields[phi_idx]->GetPhys(), first_deriv);
  m_fields[ne_idx]->PhysDeriv(2, first_deriv, HWterm_phi);
  // Compute m_HW_coeff*(d^2 n/dz^2 - d^2 Phi/dz^2)
  Array<OneD, NekDouble> HWterm_all(nPts);
  Vmath::Vsub(nPts, HWterm_ne, 1, HWterm_phi, 1, HWterm_all, 1);
  Vmath::Smul(nPts, m_HW_coeff, HWterm_all, 1, HWterm_all, 1);
  // Add result to RHS of ne and w eqns
  Vmath::Vadd(nPts, outarray[ne_idx], 1, HWterm_all, 1, outarray[ne_idx], 1);
  Vmath::Vadd(nPts, outarray[w_idx], 1, HWterm_all, 1, outarray[w_idx], 1);
}

void HWSystem::LoadParams() {
  H3LAPDSystem::LoadParams();

  // Hasegawa-Waketani coefficient (= omega/nu)
  m_session->LoadParameter("HW_coeff", m_HW_coeff, 1e4);

  // Helmholtz solve coefficient
  m_session->LoadParameter("d22", m_d22, 0);

  double tmp;
  m_session->LoadParameter("use_varcoeffs", tmp, 1);
  m_use_varcoeffs = (tmp > 0);
}

/**
 * Override SolvePhi in order to set rhs = -w, experiment with variable
 * coefficients and different values of d22
 */
void HWSystem::SolvePhi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray) {
  // Field indices
  int nPts = GetNpoints();
  int phi_idx = m_field_to_index.get_idx("phi");
  int w_idx = m_field_to_index.get_idx("w");

  // Set rhs = -w (was w * B^2 / m_d in full H3LAPD implementation)
  Array<OneD, NekDouble> rhs(nPts, 0.0);
  Vmath::Smul(nPts, -1.0, inarray[w_idx], 1, rhs, 1);

  // Set \lambda = zero (Helmholtz => Poisson)
  StdRegions::ConstFactorMap factors;
  factors[StdRegions::eFactorLambda] = 0.0;

  StdRegions::VarCoeffMap varcoeffs;
  varcoeffs[StdRegions::eVarCoeffD00] = Array<OneD, NekDouble>(nPts, 1.0);
  varcoeffs[StdRegions::eVarCoeffD01] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD02] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD11] = Array<OneD, NekDouble>(nPts, 1.0);
  varcoeffs[StdRegions::eVarCoeffD12] = Array<OneD, NekDouble>(nPts, 0.0);
  varcoeffs[StdRegions::eVarCoeffD22] = Array<OneD, NekDouble>(nPts, m_d22);

  // Solve, then backwards transform the output coefficients to update physical
  // values of Phi
  if (m_use_varcoeffs) {
    m_fields[phi_idx]->HelmSolve(rhs, m_fields[phi_idx]->UpdateCoeffs(),
                                 factors, varcoeffs);
  } else {
    m_fields[phi_idx]->HelmSolve(rhs, m_fields[phi_idx]->UpdateCoeffs(),
                                 factors);
  }
  m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(),
                              m_fields[phi_idx]->UpdatePhys());
}

} // namespace Nektar

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
// Description: 2D Hasegawa-Waketani equation system as an intermediate step
// towards the full H3-LAPD problem.  Implemented by Ed Threlfall in August 2023
// after realizing he didn't know how to do numerical flux terms in 3D
// Hasegawa-Wakatani. Parameter choices are same as in Nektar-Driftwave 2D
// proxyapp. Evolves ne, w, phi only, no momenta, no ions
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

#include "HWSystem.hpp"

namespace Nektar {
std::string HWSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "HWLAPD", HWSystem::create,
        "(2D) Hasegawa-Waketani equation system as an intermediate step "
        "towards the full H3-LAPD problem");

HWSystem::HWSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      H3LAPDSystem(pSession, pGraph) {
  m_required_flds = {"ne", "w", "phi", "ne_src"};
  m_int_fld_names = {"ne", "w"};
}

void HWSystem::AddDiffTerms(
    std::vector<std::string> field_names,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    std::vector<std::string> eqn_labels) {

  // Default is to add result of diffusing field f to the RHS of df/dt equation
  if (eqn_labels.empty()) {
    eqn_labels = std::vector(field_names);
  } else {
    ASSERTL1(field_names.size() == eqn_labels.size(),
             "AddDiffTerms: Number of quantities being diffused must match the "
             "number of equation labels.");
  }

  int nfields = field_names.size();
  int npts = inarray[0].size();

  // Make temporary copies of target fields, inarray vals and initialise a
  // temporary output array
  Array<OneD, MultiRegions::ExpListSharedPtr> tmp_fields(nfields);
  Array<OneD, Array<OneD, NekDouble>> tmp_inarray(nfields);
  Array<OneD, Array<OneD, NekDouble>> tmp_outarray(nfields);
  for (auto ii = 0; ii < nfields; ii++) {
    int idx = m_field_to_index.get_idx(field_names[ii]);
    tmp_fields[ii] = m_fields[idx];
    tmp_inarray[ii] = Array<OneD, NekDouble>(npts);
    Vmath::Vcopy(npts, inarray[idx], 1, tmp_inarray[ii], 1);
    tmp_outarray[ii] = Array<OneD, NekDouble>(outarray[idx].size());
  }
  // Compute diffusion terms; result is returned in temporary output array
  m_diffusion->Diffuse(tmp_fields.size(), tmp_fields, tmp_inarray,
                       tmp_outarray);

  // Add temporary output array to the appropriate indices of outarray
  for (auto ii = 0; ii < nfields; ii++) {
    int idx = m_field_to_index.get_idx(eqn_labels[ii]);
    Vmath::Vadd(outarray[idx].size(), outarray[idx], 1, tmp_outarray[ii], 1,
                outarray[idx], 1);
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

  // Calculate electric field from Phi, as well as corresponding drift velocity
  CalcEAndAdvVels(inarray);

  // Get field indices
  int nPts = GetNpoints();
  int ne_idx = m_field_to_index.get_idx("ne");
  int phi_idx = m_field_to_index.get_idx("phi");
  int w_idx = m_field_to_index.get_idx("w");

  // Advect ne and w (m_vAdvElec === m_vExB for HW)
  AddAdvTerms({"ne"}, m_advElec, m_vAdvElec, inarray, outarray, time);
  AddAdvTerms({"w"}, m_advVort, m_vExB, inarray, outarray, time);

  // Add \alpha*(\phi-n_e) to RHS
  Array<OneD, NekDouble> HWterm_2D_alpha(nPts);
  Vmath::Vsub(nPts, m_fields[phi_idx]->GetPhys(), 1,
              m_fields[ne_idx]->GetPhys(), 1, HWterm_2D_alpha, 1);
  Vmath::Smul(nPts, m_alpha, HWterm_2D_alpha, 1, HWterm_2D_alpha, 1);
  Vmath::Vadd(nPts, outarray[w_idx], 1, HWterm_2D_alpha, 1, outarray[w_idx], 1);
  Vmath::Vadd(nPts, outarray[ne_idx], 1, HWterm_2D_alpha, 1, outarray[ne_idx],
              1);

  // Add \kappa*\dpartial\phi/\dpartial y to RHS
  Array<OneD, NekDouble> HWterm_2D_kappa(nPts);
  m_fields[phi_idx]->PhysDeriv(1, m_fields[phi_idx]->GetPhys(),
                               HWterm_2D_kappa);
  Vmath::Smul(nPts, m_kappa, HWterm_2D_kappa, 1, HWterm_2D_kappa, 1);
  Vmath::Vsub(nPts, outarray[ne_idx], 1, HWterm_2D_kappa, 1, outarray[ne_idx],
              1);

  // Add diffusion terms
  AddDiffTerms({"ne", "w"}, inarray, outarray);

  // Add particle sources
  AddParticleSources({"ne"}, outarray);
}

/**
 * @brief Return the flux vector for the unsteady diffusion problem.
 */
void HWSystem::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor) {
  boost::ignore_unused(inarray);

  unsigned int nDim = qfield.size();
  unsigned int nDiffFields = qfield[0].size();
  unsigned int nPts = qfield[0][0].size();

  for (unsigned int j = 0; j < nDim; ++j) {
    for (unsigned int i = 0; i < nDiffFields; ++i) {
      Vmath::Smul(nPts, m_epsilon, qfield[j][i], 1, viscousTensor[j][i], 1);
    }
  }
}

// Set Phi solve RHS = w
void HWSystem::GetPhiSolveRHS(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, NekDouble> &rhs) {
  int nPts = GetNpoints();
  int w_idx = m_field_to_index.get_idx("w");
  // OP: Orig version has RHS=w, presumably it ought be alpha*w? :
  // Vmath::Smul(nPts, m_alpha, inarray[w_idx], 1, rhs, 1);
  Vmath::Vcopy(nPts, inarray[w_idx], 1, rhs, 1);
}

void HWSystem::LoadParams() {
  H3LAPDSystem::LoadParams();

  // alpha
  m_session->LoadParameter("HW_alpha", m_alpha, 2);

  // kappa
  m_session->LoadParameter("HW_kappa", m_kappa, 1);

  // diffusion type
  m_session->LoadSolverInfo("DiffusionType", m_difftype, "LDG");

  // epsilon
  m_session->LoadParameter("diffusion_coeff", m_epsilon, 0);
}

/**
 * @brief Initialization for HWSystem class.
 */
void HWSystem::v_InitObject(bool DeclareField) {
  H3LAPDSystem::v_InitObject(DeclareField);

  // Setup diffusion object
  std::string diffName;
  m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
  m_diffusion =
      SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
  m_diffusion->SetFluxVector(&HWSystem::GetFluxVectorDiff, this);
  m_diffusion->InitObject(m_session, m_fields);
}

} // namespace Nektar

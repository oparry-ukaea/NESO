///////////////////////////////////////////////////////////////////////////////
//
// File Outflow1DSystem.cpp
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

#include "Outflow1DSystem.h"

namespace Nektar {
std::string Outflow1DSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "Outflow1D", Outflow1DSystem::create, "Simple 1D outflow problem");

Outflow1DSystem::Outflow1DSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      H3LAPDSystem(pSession, pGraph) {
  m_required_flds = {"ne", "Ge"};
}

void Outflow1DSystem::ExplicitTimeInt(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time) {

  // Check inarray for NaNs
  for (auto &var : {"ne", "Ge"}) {
    auto fidx = m_field_to_index.get_idx(var);
    for (auto ii = 0; ii < inarray[fidx].size(); ii++) {
      if (!std::isfinite(inarray[fidx][ii])) {
        std::cout << "NaN in field " << var << ", aborting." << std::endl;
        exit(1);
      }
    }
  }

  ZeroOutArray(outarray);

  // Compute ve from Ge

  // Get field indices
  int nPts = GetNpoints();
  int ne_idx = m_field_to_index.get_idx("ne");
  int Ge_idx = m_field_to_index.get_idx("Ge");

  // v_par,e = Ge / max(ne,n_floor) / me
  int Ge_idx = m_field_to_index.get_idx("Ge");
  for (auto idim = 0; idim < m_graph->GetSpaceDimension(); idim++) {
    for (auto ii = 0; ii < nPts; ii++) {
      m_advElec[idim][ii] =
          inarray[Ge_idx][ii] /
          std::max(inarray[ne_idx][ii], m_nRef * m_n_floor_fac);
    }
  }

  AddAdvTerms({"ne"}, m_advElec, m_vExB, inarray, outarray, time);

  // // Advect both ne and w using ExB velocity
  // AddAdvTerms({"ne"}, m_advElec, m_vExB, inarray, outarray, time);
  // AddAdvTerms({"w"}, m_advVort, m_vExB, inarray, outarray, time);

  // // Add \alpha*(\phi-n_e) to RHS
  // Array<OneD, NekDouble> HWterm_2D_alpha(nPts);
  // Vmath::Vsub(nPts, m_fields[phi_idx]->GetPhys(), 1,
  //             m_fields[ne_idx]->GetPhys(), 1, HWterm_2D_alpha, 1);
  // Vmath::Smul(nPts, m_alpha, HWterm_2D_alpha, 1, HWterm_2D_alpha, 1);
  // Vmath::Vadd(nPts, outarray[w_idx], 1, HWterm_2D_alpha, 1, outarray[w_idx],
  // 1); Vmath::Vadd(nPts, outarray[ne_idx], 1, HWterm_2D_alpha, 1,
  // outarray[ne_idx],
  //             1);

  // // Add \kappa*\dpartial\phi/\dpartial y to RHS
  // Array<OneD, NekDouble> HWterm_2D_kappa(nPts);
  // m_fields[phi_idx]->PhysDeriv(1, m_fields[phi_idx]->GetPhys(),
  //                              HWterm_2D_kappa);
  // Vmath::Vsub(nPts, outarray[ne_idx], 1, HWterm_2D_kappa, 1,
  // outarray[ne_idx],
  //             1);
}

void Outflow1DSystem::LoadParams() {
  H3LAPDSystem::LoadParams();

  // // c
  // m_session->LoadParameter("c", m_c);

  // // R
  // m_session->LoadParameter("R", m_R);
}

} // namespace Nektar

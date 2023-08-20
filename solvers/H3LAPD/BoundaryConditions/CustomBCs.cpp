///////////////////////////////////////////////////////////////////////////////
//
// File: CustomBCs.cpp
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
// Description: Abstract base class for custom boundary conditions.
//
///////////////////////////////////////////////////////////////////////////////

#include "CustomBCs.h"

using namespace std;

namespace Nektar {
CustomBCsFactory &GetCustomBCsFactory() {
  static CustomBCsFactory instance;
  return instance;
}

CustomBCs::CustomBCs(const LibUtilities::SessionReaderSharedPtr &session,
                     const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                     const Array<OneD, Array<OneD, NekDouble>> &trace_normals,
                     const int field_idx, const int space_dim,
                     const int bc_region, const int cnt,
                     SpatialDomains::BoundaryConditionShPtr cnd)
    : m_session(session), m_fields(fields), m_traceNormals(trace_normals),
      m_spacedim(space_dim), m_bcRegion(bc_region), m_offset(cnt), m_cnd(cnd),
      m_field_to_index(session->GetVariables()) {
  ASSERTL0(space_dim == 3, "CustomBCs only set up for space_dim=3");
  SpatialDomains::BoundaryConditionType BCtype =
      m_cnd->GetBoundaryConditionType();
  ASSERTL0(BCtype == SpatialDomains::eDirichlet ||
               BCtype == SpatialDomains::eNeumann,
           "CustomBCs must be either Dirichlet or Neumann type");
}

void CustomBCs::evaluate_expression(LocalRegions::ExpansionSharedPtr explist,
                                    Array<OneD, NekDouble> result) {

  int npts = explist->GetTotPoints();
  // Coords of quad points in this explist
  Array<OneD, NekDouble> tmp_x(npts), tmp_y(npts), tmp_z(npts);
  explist->GetCoords(tmp_x, tmp_y, tmp_y);

  switch (m_cnd->GetBoundaryConditionType()) {
  case SpatialDomains::eDirichlet: {
    auto dcond =
        std::dynamic_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(
            m_cnd);
    dcond->m_dirichletCondition.Evaluate(tmp_x, tmp_y, tmp_y, result);
  } break;
  case SpatialDomains::eNeumann: {
    auto ncond =
        std::dynamic_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(
            m_cnd);
    ncond->m_neumannCondition.Evaluate(tmp_x, tmp_y, tmp_y, result);
  } break;
  default: {
  } break;
  }
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void CustomBCs::Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                      Array<OneD, Array<OneD, NekDouble>> &physarray,
                      const NekDouble &time) {
  v_Apply(Fwd, physarray, time);
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value, like Dirichlet,
 * bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions
 * is the target value like WallViscousBC weight should be 0.5.
 */
void CustomBCs::v_ApplyBwdWeight() {
  for (int i = 0; i < m_fields.size(); ++i) {
    m_fields[i]->SetBndCondBwdWeight(m_bcRegion, m_weight);
  }
}

} // namespace Nektar

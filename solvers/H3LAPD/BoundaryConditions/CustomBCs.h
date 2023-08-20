///////////////////////////////////////////////////////////////////////////////
//
// File: CustomBCs.h
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

#ifndef H3LAPD_CUSTOMBCS_H
#define H3LAPD_CUSTOMBCS_H

#include "nektar_interface/utilities.hpp"

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Conditions.h>
#include <string>

namespace Nektar {
//  Forward declaration
class CustomBCs;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<CustomBCs> CustomBCsSharedPtr;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<
    std::string, CustomBCs, const LibUtilities::SessionReaderSharedPtr &,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &,
    const Array<OneD, Array<OneD, NekDouble>> &, const int, const int,
    const int, const int, SpatialDomains::BoundaryConditionShPtr>
    CustomBCsFactory;

/// Declaration of the boundary condition factory singleton
CustomBCsFactory &GetCustomBCsFactory();

/**
 * @class CustomBCs.
 */
class CustomBCs {
public:
  virtual ~CustomBCs() {}

  /// Apply the boundary condition
  void Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
             Array<OneD, Array<OneD, NekDouble>> &physarray,
             const NekDouble &time = 0);

  /// Apply the Weight of boundary condition
  void ApplyBwdWeight() { v_ApplyBwdWeight(); }

protected:
  /// Session reader
  LibUtilities::SessionReaderSharedPtr m_session;
  /// Array of fields
  Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
  /// Trace normals
  Array<OneD, Array<OneD, NekDouble>> m_traceNormals;
  // Field name => index mapper
  NESO::NektarFieldIndexMap m_field_to_index;
  /// Space dimension
  int m_spacedim;
  /// Id of the boundary region
  int m_bcRegion;
  /// Offset
  int m_offset;
  /// Boundary weight
  NekDouble m_weight;

  // Nektar BC ptr
  SpatialDomains::BoundaryConditionShPtr m_cnd;

  /// Constructor
  CustomBCs(const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
            const int field_idx, const int pSpaceDim, const int bcRegion,
            const int cnt, SpatialDomains::BoundaryConditionShPtr cnd);

  void evaluate_expression(LocalRegions::ExpansionSharedPtr explist,
                           Array<OneD, NekDouble> result);

  virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                       Array<OneD, Array<OneD, NekDouble>> &physarray,
                       const NekDouble &time) = 0;

  virtual void v_ApplyBwdWeight();
};
} // namespace Nektar

#endif

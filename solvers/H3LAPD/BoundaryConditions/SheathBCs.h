///////////////////////////////////////////////////////////////////////////////
//
// File: SheathBCs.h
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
// Description: Sheath boundary conditions class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef H3LAPD_SHEATHBCS_H
#define H3LAPD_SHEATHBCS_H

#include "CustomBCs.h"

namespace Nektar {

class SheathBCs : public CustomBCs {
public:
  friend class MemoryManager<SheathBCs>;

  /// Creates an instance of this class
  static CustomBCsSharedPtr
  create(const LibUtilities::SessionReaderSharedPtr &pSession,
         const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
         const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
         const int field_idx, const int pSpaceDim, const int bcRegion,
         const int cnt, SpatialDomains::BoundaryConditionShPtr cnd) {
    CustomBCsSharedPtr p = MemoryManager<SheathBCs>::AllocateSharedPtr(
        pSession, pFields, pTraceNormals, field_idx, pSpaceDim, bcRegion, cnt,
        cnd);
    return p;
  }

  /// Name of the class
  static std::string className;

protected:
  virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                       Array<OneD, Array<OneD, NekDouble>> &physarray,
                       const NekDouble &time);

private:
  SheathBCs(const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
            const int field_idx, const int pSpaceDim, const int bcRegion,
            const int cnt, SpatialDomains::BoundaryConditionShPtr cnd);

  virtual ~SheathBCs(void){};

  int m_dens_idx;
  int m_mom_idx;
};

} // namespace Nektar

#endif

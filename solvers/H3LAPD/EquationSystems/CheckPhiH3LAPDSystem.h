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
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CheckPhiH3LAPDSystem_H
#define CheckPhiH3LAPDSystem_H

#include "nektar_interface/utilities.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include "H3LAPDSystem.h"

namespace Nektar {

class CheckPhiH3LAPDSystem : virtual public H3LAPDSystem {
public:
  friend class MemoryManager<CheckPhiH3LAPDSystem>;

  /// Name of class.
  static std::string className;

  /// Creates an instance of this class.
  static SolverUtils::EquationSystemSharedPtr
  create(const LibUtilities::SessionReaderSharedPtr &pSession,
         const SpatialDomains::MeshGraphSharedPtr &pGraph) {
    SolverUtils::EquationSystemSharedPtr p =
        MemoryManager<CheckPhiH3LAPDSystem>::AllocateSharedPtr(pSession,
                                                               pGraph);
    p->InitObject();
    return p;
  }

protected:
  CheckPhiH3LAPDSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr &pGraph);

  void ExplicitTimeInt(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time) override;

  void GetPhiSolveRHS(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                      Array<OneD, NekDouble> &rhs) override;
};

} // namespace Nektar
#endif

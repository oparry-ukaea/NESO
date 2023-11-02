///////////////////////////////////////////////////////////////////////////////
//
// File: CustomUpwindSolver.h
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
// Description: Upwind Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef H3LAPD_CUSTOM_UPWIND_SOLVER_H
#define H3LAPD_CUSTOM_UPWIND_SOLVER_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;

namespace NESO::Solvers::H3LAPD {

typedef Nektar::Array<Nektar::OneD, Nektar::NekDouble> Nek1DArr;
typedef Nektar::Array<Nektar::OneD, Nek1DArr> Nek2DArr;
typedef Nektar::Array<Nektar::OneD, const Nek1DArr> Nek2DArrConstInner;

/**
 * @brief Subclass Nektar's Riemann solver in order to set up non-standard
 * fluxes.
 */
class CustomUpwindSolver : public SU::RiemannSolver {
public:
  SOLVER_UTILS_EXPORT static SU::RiemannSolverSharedPtr
  create(const LU::SessionReaderSharedPtr &pSession) {
    return SU::RiemannSolverSharedPtr(new CustomUpwindSolver(pSession));
  }

  static std::string solver_name;

protected:
  CustomUpwindSolver(const LU::SessionReaderSharedPtr &pSession);

  virtual void v_Solve(const int nDim, const Nek2DArrConstInner &Fwd,
                       const Nek2DArrConstInner &Bwd,
                       Nek2DArr &flux) override final;
};
} // namespace NESO::Solvers::H3LAPD

#endif

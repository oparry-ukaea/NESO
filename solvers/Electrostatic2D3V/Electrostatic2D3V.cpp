///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSolver.cpp
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
// Description: Compressible Flow Equations framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Timer.h>

#include <SolverUtils/Driver.h>
#include <SolverUtils/EquationSystem.h>
#include <SpatialDomains/MeshGraphIO.h>

#include "ElectrostaticTwoStream2D3V.hpp"

#include <memory>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

#ifndef FIELD_TYPE
#define FIELD_TYPE ContField
#endif

int main(int argc, char *argv[]) {
  LibUtilities::SessionReaderSharedPtr session;
  SpatialDomains::MeshGraphSharedPtr graph;
  string vDriverModule;
  DriverSharedPtr drv;

  try {
    // Create session reader.
    session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    // Create MeshGraph.
    graph = SpatialDomains::MeshGraphIO::Read(session);

    // Create driver
    session->LoadSolverInfo("Driver", vDriverModule, "Standard");
    drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

    auto electrostatic_two_stream_2d3v =
        std::make_shared<ElectrostaticTwoStream2D3V<FIELD_TYPE>>(session, graph,
                                                                 drv);
    electrostatic_two_stream_2d3v->run();

    // Print out timings if verbose
    if (session->DefinesCmdLineArgument("verbose")) {
      int iolevel;

      session->LoadParameter("IO_Timer_Level", iolevel, 1);

      LibUtilities::Timer::PrintElapsedRegions(session->GetComm(), std::cout,
                                               iolevel);
    }

    electrostatic_two_stream_2d3v->finalise();
    // Finalise communications
    session->Finalise();
  } catch (const std::runtime_error &) {
    return 1;
  } catch (const std::string &eStr) {
    cout << "Error: " << eStr << endl;
  }

  return 0;
}

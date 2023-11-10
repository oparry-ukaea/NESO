///////////////////////////////////////////////////////////////////////////////
//
// File: SimpleSOL.cpp
//
//
// Description: Driver for SimpleSOL.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Driver.h>
#include <limits>

#include "SimpleSOL.h"

using namespace Nektar;
using namespace Nektar::SolverUtils;

namespace NESO {
namespace Solvers {

//======================= Helper functions ====================================
/**
 * Computes min,max x,y for a vector of point PointGeom objects
 */
void compute_minmax_xy(
    const std::vector<SpatialDomains::PointGeomSharedPtr> &verts, double &xmin,
    double &xmax, double &ymin, double &ymax) {
  for (auto vert : verts) {
    Array<OneD, NekDouble> pcoords;
    vert->GetCoords(pcoords);
    xmin = std::min(xmin, pcoords[0]);
    xmax = std::max(xmax, pcoords[0]);
    ymin = std::min(ymin, pcoords[1]);
    ymax = std::max(ymax, pcoords[1]);
  }
}

/**
 * Extracts the size of the longest geometry vector in a composite map
 */
std::size_t max_composite_size(const SpatialDomains::CompositeMap &composites,
                               LibUtilities::ShapeType type) {
  return std::transform_reduce(
      composites.cbegin(), composites.cend(), std::size_t{0},
      [](const auto &a, const auto &b) { return std::max(a, b); },
      [type](auto c) {
        return (c.second->m_geomVec[0]->GetShapeType() == type)
                   ? c.second->m_geomVec.size()
                   : 0;
      });
}

void extract_mesh_props(LibUtilities::SessionReaderSharedPtr session,
                        SpatialDomains::MeshGraphSharedPtr graph,
                        double &mesh_tan_theta, double &mesh_max_dx) {
  SpatialDomains::CompositeMap composites = graph->GetComposites();
  std::vector<SpatialDomains::PointGeomSharedPtr> edge_verts;

  switch (graph->GetMeshDimension()) {
  case 1:
    break;
  case 2:
    std::size_t long_dim_nsegs =
        max_composite_size(composites, LibUtilities::ShapeType::eSegment);
    for (auto c : composites) {

      SpatialDomains::GeometrySharedPtr first = c.second->m_geomVec[0];
      if (first->GetShapeType() == LibUtilities::ShapeType::eSegment) {
        std::size_t nsegs = c.second->m_geomVec.size();
        // std::cout << "Comp " << c.first << " (len = " << nsegs << ")"
        //           << std::endl;
        SpatialDomains::GeometrySharedPtr last = c.second->m_geomVec[nsegs - 1];
        // Compute dx,dy using vertices belonging to end segments
        edge_verts = {first->GetVertex(0), first->GetVertex(1),
                      last->GetVertex(0), last->GetVertex(1)};
        double xmin(std::numeric_limits<double>::max()), ymin = xmin;
        double xmax(std::numeric_limits<double>::min()), ymax = xmax;
        compute_minmax_xy(edge_verts, xmin, xmax, ymin, ymax);
        session->GetComm()->AllReduce(xmin, LibUtilities::ReduceMin);
        if (nsegs == long_dim_nsegs) {
        } else {
        }
      }
    }

    break;
  }
}

/**
 * Check that the parameters theta and s_max are consistent with the mesh size
 * and orientation.
 */
int check_mesh(LibUtilities::SessionReaderSharedPtr session,
               SpatialDomains::MeshGraphSharedPtr graph,
               const double tol = 1e-6) {
  // Get params from session
  double theta, smax;
  session->LoadParameter("theta", theta, 0.0);
  session->LoadParameter("srcs_smax", smax, 110.0);
  double tan_theta_config = std::tan(theta);

  double tan_theta_mesh, mesh_max_dx;
  extract_mesh_props(session, graph, tan_theta_mesh, mesh_max_dx);

  // Check size of long dimension
  double long_dim_size = smax * std::cos(theta);
  double x_diff = std::abs(mesh_max_dx - long_dim_size);
  if (x_diff > tol) {
    std::cerr << "Longest dimension of mesh has size=" << dx
              << "; but config's theta,s_max combination requires size "
              << long_dim_size << " (diff = " << x_diff << ") Aborting..."
              << std::endl;
    return 1;
  }

  // Check orientation along both dimensions
  double tan_theta_diff = std::abs(tan_theta_mesh - tan_theta_config);
  if (tan_theta_diff > tol) {
    std::cerr << "Mesh implies tan(theta)=" << tan_theta_mesh
              << "; but config sets tan(theta)=" << tan_theta_config
              << " (diff = " << tan_theta_diff << ") Aborting..." << std::endl;
    return 2;
  }
  return 0;
}
//=============================================================================

int run_SimpleSOL(int argc, char *argv[]) {
  try {
    // Create session reader.
    auto session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    // Read the mesh and create a MeshGraph object.
    auto graph = SpatialDomains::MeshGraph::Read(session);

    // Check mesh, config are consistent
    if (check_mesh(session, graph)) {
      session->Finalise();
      return 4;
    }

    // Create driver.
    std::string driverName;
    session->LoadSolverInfo("Driver", driverName, "Standard");
    auto drv = GetDriverFactory().CreateInstance(driverName, session, graph);

    // Execute driver
    drv->Execute();

    // Finalise session
    session->Finalise();
  } catch (const std::runtime_error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (const std::string &eStr) {
    std::cerr << "Error: " << eStr << std::endl;
    return 2;
  }

  return 0;
}

} // namespace Solvers
} // namespace NESO
#ifndef __COMPOSITE_INTERSECTION_H_
#define __COMPOSITE_INTERSECTION_H_

#include <SpatialDomains/MeshGraph.h>
using namespace Nektar;

#include <neso_particles.hpp>
using namespace NESO::Particles;

#include <nektar_interface/geometry_transport/packed_geom_2d.hpp>
#include <nektar_interface/particle_mesh_interface.hpp>

#include "composite_transport.hpp"

#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace NESO::CompositeInteraction {

/**
 *  High-level class to detect and compute the intersection of a particle
 *  trajectory and a Nektar++ composite.
 */
class CompositeIntersection {
protected:
  const int ndim;
  ParticleMeshInterfaceSharedPtr particle_mesh_interface;
  std::unique_ptr<CompositeTransport> composite_transport;

public:
  /// Disable (implicit) copies.
  CompositeIntersection(const CompositeIntersection &st) = delete;
  /// Disable (implicit) copies.
  CompositeIntersection &operator=(CompositeIntersection const &a) = delete;

  /// The NESO::Particles Sym<REAL> used to store the previous particle
  /// position.
  const static inline Sym<REAL> previous_position_sym =
      Sym<REAL>("NESO_COMP_INT_PREV_POS");

  /// The composite indices for which the class detects intersections with.
  const std::vector<int> composite_indices;

  /**
   * TODO
   */
  inline void free() { this->composite_transport->free(); }

  /**
   *  TODO
   */
  CompositeIntersection(ParticleMeshInterfaceSharedPtr particle_mesh_interface,
                        std::vector<int> &composite_indices)
      : particle_mesh_interface(particle_mesh_interface),
        ndim(particle_mesh_interface->graph->GetMeshDimension()),
        composite_indices(composite_indices) {

    this->composite_transport = std::make_unique<CompositeTransport>(
        particle_mesh_interface, composite_indices);
  }

  /**
   *  Method to store the current particle positions before an integration step.
   *
   *  @param particle_group Particles to store current positions of.
   */
  inline void pre_integration(ParticleGroupSharedPtr particle_group) {
    const auto position_dat = particle_group->position_dat;
    const int ndim = position_dat->ncomp;
    NESOASSERT(ndim == this->ndim,
               "missmatch between particle ndim and class ndim");
    const auto sycl_target = particle_group->sycl_target;
    // If the previous position dat does not already exist create it here
    if (!particle_group->contains_dat(previous_position_sym)) {
      particle_group->add_particle_dat(
          ParticleDat(sycl_target, ParticleProp(previous_position_sym, ndim),
                      particle_group->domain->mesh->get_cell_count()));
    }

    // copy the current position onto the previous position
    auto pl_iter_range = position_dat->get_particle_loop_iter_range();
    auto pl_stride = position_dat->get_particle_loop_cell_stride();
    auto pl_npart_cell = position_dat->get_particle_loop_npart_cell();
    const auto k_P = position_dat->cell_dat.device_ptr();
    auto k_PP =
        particle_group->get_dat(previous_position_sym)->cell_dat.device_ptr();
    sycl_target->queue
        .submit([&](sycl::handler &cgh) {
          cgh.parallel_for<>(
              sycl::range<1>(pl_iter_range), [=](sycl::id<1> idx) {
                NESO_PARTICLES_KERNEL_START
                const INT cellx = NESO_PARTICLES_KERNEL_CELL;
                const INT layerx = NESO_PARTICLES_KERNEL_LAYER;
                for (int dimx = 0; dimx < ndim; dimx++) {
                  k_PP[cellx][dimx][layerx] = k_P[cellx][dimx][layerx];
                }
                NESO_PARTICLES_KERNEL_END
              });
        })
        .wait_and_throw();
  }
};

} // namespace NESO::CompositeInteraction

#endif

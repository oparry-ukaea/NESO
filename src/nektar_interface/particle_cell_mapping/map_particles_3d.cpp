#include <nektar_interface/particle_cell_mapping/map_particles_3d.hpp>

namespace NESO {

MapParticles3D::MapParticles3D(
    SYCLTargetSharedPtr sycl_target,
    ParticleMeshInterfaceSharedPtr particle_mesh_interface,
    ParameterStoreSharedPtr config)
    : sycl_target(sycl_target),
      particle_mesh_interface(particle_mesh_interface) {

  this->map_particles_common =
      std::make_unique<MapParticlesCommon>(sycl_target);

  this->map_particles_3d_regular = nullptr;
  this->map_particles_host = nullptr;
  std::get<0>(this->map_particles_3d_deformed_linear) = nullptr;
  std::get<1>(this->map_particles_3d_deformed_linear) = nullptr;
  std::get<2>(this->map_particles_3d_deformed_linear) = nullptr;
  std::get<3>(this->map_particles_3d_deformed_linear) = nullptr;
  this->map_particles_3d_deformed_non_linear = nullptr;

  GeometryContainer3D geometry_container_3d;
  const bool all_generic_newton = static_cast<bool>(
      config->get<INT>("MapParticles3D/all_generic_newton", 0));

  if (!all_generic_newton) {
    // Get the geometry objects in their different categories.
    assemble_geometry_container_3d(particle_mesh_interface->graph,
                                   particle_mesh_interface->remote_geoms_3d,
                                   geometry_container_3d);
  } else {
    // If we are putting all geoms through the generic 3D mapper then only
    // populate the deformed_non_linear container.
    {
      std::map<int, std::shared_ptr<Nektar::SpatialDomains::Geometry3D>> geoms;
      get_all_elements_3d(particle_mesh_interface->graph, geoms);
      for (auto gx : geoms) {
        std::pair<int, std::shared_ptr<Nektar::SpatialDomains::Geometry3D>>
            tmp = {gx.first, gx.second};
        geometry_container_3d.deformed_non_linear.push_back(tmp);
      }
    }
    {
      for (auto gx : particle_mesh_interface->remote_geoms_3d) {
        geometry_container_3d.deformed_non_linear.push_back(gx);
      }
    }
  }

  // Create a mapper for 3D regular geometry objects
  if (geometry_container_3d.regular.size()) {
    this->map_particles_3d_regular = std::make_unique<MapParticles3DRegular>(
        sycl_target, particle_mesh_interface, config);
  }

  // Create mappers for the deformed geometry objects with linear faces
  if (geometry_container_3d.deformed_linear.tet.size()) {
    std::get<0>(this->map_particles_3d_deformed_linear) = std::make_unique<
        Newton::MapParticlesNewton<Newton::MappingTetLinear3D>>(
        Newton::MappingTetLinear3D{}, this->sycl_target,
        geometry_container_3d.deformed_linear.tet.local,
        geometry_container_3d.deformed_linear.tet.remote, config);
  }
  if (geometry_container_3d.deformed_linear.prism.size()) {
    std::get<1>(this->map_particles_3d_deformed_linear) = std::make_unique<
        Newton::MapParticlesNewton<Newton::MappingPrismLinear3D>>(
        Newton::MappingPrismLinear3D{}, this->sycl_target,
        geometry_container_3d.deformed_linear.prism.local,
        geometry_container_3d.deformed_linear.prism.remote, config);
  }
  if (geometry_container_3d.deformed_linear.hex.size()) {
    std::get<2>(this->map_particles_3d_deformed_linear) = std::make_unique<
        Newton::MapParticlesNewton<Newton::MappingHexLinear3D>>(
        Newton::MappingHexLinear3D{}, this->sycl_target,
        geometry_container_3d.deformed_linear.hex.local,
        geometry_container_3d.deformed_linear.hex.remote, config);
  }
  if (geometry_container_3d.deformed_linear.pyr.size()) {
    std::get<3>(this->map_particles_3d_deformed_linear) = std::make_unique<
        Newton::MapParticlesNewton<Newton::MappingPyrLinear3D>>(
        Newton::MappingPyrLinear3D{}, this->sycl_target,
        geometry_container_3d.deformed_linear.pyr.local,
        geometry_container_3d.deformed_linear.pyr.remote, config);
  }
  if (geometry_container_3d.deformed_non_linear.size()) {

    std::map<int, std::shared_ptr<Geometry3D>> local;
    std::vector<std::shared_ptr<RemoteGeom3D>> remote;
    remote.reserve(
        geometry_container_3d.deformed_non_linear.tet.remote.size() +
        geometry_container_3d.deformed_non_linear.pyr.remote.size() +
        geometry_container_3d.deformed_non_linear.prism.remote.size() +
        geometry_container_3d.deformed_non_linear.hex.remote.size());

    auto lambda_push = [&](auto &lr) -> void {
      for (auto &lx : lr.local) {
        local[lx.first] = lx.second;
      }
      for (auto &rx : lr.remote) {
        remote.push_back(rx);
      }
    };
    lambda_push(geometry_container_3d.deformed_non_linear.tet);
    lambda_push(geometry_container_3d.deformed_non_linear.pyr);
    lambda_push(geometry_container_3d.deformed_non_linear.prism);
    lambda_push(geometry_container_3d.deformed_non_linear.hex);

    this->map_particles_3d_deformed_non_linear =
        std::make_unique<Newton::MapParticlesNewton<Newton::MappingGeneric3D>>(
            Newton::MappingGeneric3D{}, this->sycl_target, local, remote,
            config);
  }

  // Create a host mapper as a last resort mapping attempt.
  this->map_particles_host = std::make_unique<MapParticlesHost>(
      sycl_target, particle_mesh_interface, config);
}

void MapParticles3D::map(ParticleGroup &particle_group, const int map_cell) {

  if (this->map_particles_3d_regular) {
    // attempt to bin particles into regular geometry objects
    this->map_particles_3d_regular->map(particle_group, map_cell);
  }

  map_newton_initial(std::get<0>(this->map_particles_3d_deformed_linear),
                     particle_group, map_cell);
  map_newton_initial(std::get<1>(this->map_particles_3d_deformed_linear),
                     particle_group, map_cell);
  map_newton_initial(std::get<2>(this->map_particles_3d_deformed_linear),
                     particle_group, map_cell);
  map_newton_initial(std::get<3>(this->map_particles_3d_deformed_linear),
                     particle_group, map_cell);
  map_newton_initial(this->map_particles_3d_deformed_non_linear, particle_group,
                     map_cell);

  map_newton_final(std::get<0>(this->map_particles_3d_deformed_linear),
                   particle_group, map_cell);
  map_newton_final(std::get<1>(this->map_particles_3d_deformed_linear),
                   particle_group, map_cell);
  map_newton_final(std::get<2>(this->map_particles_3d_deformed_linear),
                   particle_group, map_cell);
  map_newton_final(std::get<3>(this->map_particles_3d_deformed_linear),
                   particle_group, map_cell);
  map_newton_final(this->map_particles_3d_deformed_non_linear, particle_group,
                   map_cell);

  if (map_cell > -1) {
    // if there are particles not yet mapped this may be an error depending on
    // which stage of NESO-Particles hybrid move we are at.
    bool particles_not_mapped =
        this->map_particles_common->check_map(particle_group, map_cell);

    if (this->map_particles_host && particles_not_mapped) {
      auto pr = ProfileRegion("MapParticles3D", "host_backup");
      this->map_particles_host->map(particle_group, -1);
      pr.end();
      this->sycl_target->profile_map.add_region(pr);
    }

    particles_not_mapped =
        this->map_particles_common->check_map(particle_group);

    if (particles_not_mapped) {
      nprint(
          "==================================================================="
          "=============");
      auto mpi_rank_dat = particle_group.mpi_rank_dat;
      auto cell_id_dat = particle_group.cell_id_dat;
      auto positions_dat = particle_group.position_dat;
      const auto cell_count = cell_id_dat->cell_dat.ncells;
      auto sym_mpi_ranks = mpi_rank_dat->sym;
      auto sym_cell = cell_id_dat->sym;
      auto sym_positions = positions_dat->sym;
      const auto ndim = particle_group.domain->mesh->get_ndim();
      nprint("rank:", this->sycl_target->comm_pair.rank_parent);
      for (int cx = 0; cx < cell_count; cx++) {
        auto P = particle_group.get_cell(sym_positions, cx);
        auto C = particle_group.get_cell(sym_cell, cx);
        auto M = particle_group.get_cell(sym_mpi_ranks, cx);
        const auto nrow = P->nrow;
        for (int rx = 0; rx < nrow; rx++) {
          if (M->at(rx, 1) < 0) {
            std::cout << std::setprecision(18);
            particle_group.print_particle(cx, rx);
          }
        }
      }
      nprint(
          "==================================================================="
          "=============");
    }

    NESOASSERT(!particles_not_mapped,
               "Failed to find cell containing one or more particles.");
  }
}

} // namespace NESO

include(CMakeDependentOption)

set(EXEC_TARGET_NAME MaxwellWave)

# create the solver library
set(MAXWELLWAVE_LIBRARY_NAME
    MaxwellWave_ObjLib
    CACHE INTERNAL "")
add_library(${MAXWELLWAVE_LIBRARY_NAME} OBJECT
            ${CMAKE_CURRENT_LIST_DIR}/EquationSystems/MaxwellWaveSystem.cpp)

target_include_directories(
  ${MAXWELLWAVE_LIBRARY_NAME}
  PRIVATE ${MPI_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/include
          ${NESO_PARTICLES_INCLUDE_PATH})

target_link_libraries(
  ${MAXWELLWAVE_LIBRARY_NAME}
  PRIVATE ${NESO_LIBRARY_NAME} ${NESO_LIBRARY_NAME} Nektar++::nektar++
          NESO-Particles::NESO-Particles)

# create the solver binary
file(GLOB MaxwellWaveSource ${CMAKE_CURRENT_LIST_DIR}/MaxwellWave2D3V.cpp)

add_executable(${EXEC_TARGET_NAME} ${MaxwellWaveSource})

target_link_libraries(
  ${EXEC_TARGET_NAME} PRIVATE ${MAXWELLWAVE_LIBRARY_NAME} Nektar++::nektar++
                              NESO-Particles::NESO-Particles)

add_link_options(${EXEC_TARGET_NAME} PRIVATE ${NESO_PARTICLES_LINK_FLAGS})

target_include_directories(
  ${EXEC_TARGET_NAME} PRIVATE ${MPI_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/include
                              ${NESO_PARTICLES_INCLUDE_PATH})

add_sycl_to_target(TARGET ${EXEC_TARGET_NAME} SOURCES ${MaxwellWaveSource})

# Install location
install(TARGETS ${EXEC_TARGET_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

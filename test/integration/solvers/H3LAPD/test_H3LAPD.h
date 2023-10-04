#ifndef H3LAPD_TESTS_H
#define H3LAPD_TESTS_H

#include <gtest/gtest.h>

#include "EquationSystems/H3LAPDSystem.hpp"
#include "H3LAPD.hpp"
#include "solver_test_utils.h"
#include "solvers/solver_callback_handler.hpp"
#include "solvers/solver_runner.hpp"

// Mass conservation tolerance
constexpr double mass_cons_tolerance = 1e-14;

struct H3LAPDMassConservationPre : public NESO::SolverCallback<H3LAPDSystem> {
  void call(H3LAPDSystem *state) {
    state->m_diag_mass_recording->compute_initial_fluid_mass();
  }
};

struct H3LAPDMassConservationPost : public NESO::SolverCallback<H3LAPDSystem> {
  std::vector<double> mass_error;
  void call(H3LAPDSystem *state) {
    auto md = state->m_diag_mass_recording;
    const double mass_particles = md->compute_particle_mass();
    const double mass_fluid = md->compute_fluid_mass();
    const double mass_total = mass_particles + mass_fluid;
    const double mass_added = md->compute_total_added_mass();
    const double correct_total = mass_added + md->get_initial_mass();
    this->mass_error.push_back(std::fabs(correct_total - mass_total) /
                               std::fabs(correct_total));
  }
};

class H3LAPDTest : public NektarSolverTest {
protected:
  void check_mass_conservation(const double &tolerance) {

    H3LAPDMassConservationPre callback_pre;
    H3LAPDMassConservationPost callback_post;

    MainFuncType runner = [&](int argc, char **argv) {
      SolverRunner solver_runner(argc, argv);
      auto equation_system = std::dynamic_pointer_cast<H3LAPDSystem>(
          solver_runner.driver->GetEqu()[0]);

      equation_system->m_solver_callback_handler.register_pre_integrate(
          callback_pre);
      equation_system->m_solver_callback_handler.register_post_integrate(
          callback_post);

      solver_runner.execute();
      solver_runner.finalise();
      return 0;
    };

    int ret_code = run(runner);
    EXPECT_EQ(ret_code, 0);
    ASSERT_THAT(callback_post.mass_error,
                testing::Each(testing::Le(mass_cons_tolerance)));
  }
};

#endif // H3LAPD_TESTS_H

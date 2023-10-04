#include <gtest/gtest.h>

#include "H3LAPD.hpp"
#include "test_H3LAPD.h"

/**
 * Tests for H3LAPD solver. Note that the test name itself is used to
 * determine the locations of the config file, mesh and initial conditions in
 * each case.
 */

TEST_F(H3LAPDTest, HWFluidOnly) {
  check_mass_conservation(mass_cons_tolerance);
}

TEST_F(H3LAPDTest, HWCoupled) { check_mass_conservation(mass_cons_tolerance); }
#include "gtest/gtest.h"
#include "gtest-mpi-listener.hpp"
#include <iostream>
#include <mpi.h>

/*
 *  If an exception is thrown try and abort MPI cleanly to prevent a deadlock.
 */
void TerminateHandler() {
  std::cout << "Exception thrown attempting to abort cleanly." << std::endl;
  MPI_Abort(MPI_COMM_WORLD, -2);
}

int main(int argc, char** argv) {
  // Filter out Google Test arguments
  #if GTEST_HAS_EXCEPTIONS
  std::set_terminate(&TerminateHandler);
#endif
  ::testing::InitGoogleTest(&argc, argv);

  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cout << "ERROR: MPI_Init != MPI_SUCCESS" << std::endl;
    return -1;
  }

  // Add object that will finalize MPI on exit; Google Test owns this pointer
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);

  // Get the event listener list.
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  // Remove default listener: the default printer and the default XML printer
  ::testing::TestEventListener *l =
        listeners.Release(listeners.default_result_printer());

  // Adds MPI listener; Google Test owns this pointer
  listeners.Append(
      new GTestMPIListener::MPIWrapperPrinter(l,
                                              MPI_COMM_WORLD)
      );
  // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
  // pass and 1 if some test fails.
  int result = RUN_ALL_TESTS();

  return 0;  // Run tests, then clean up and exit
}
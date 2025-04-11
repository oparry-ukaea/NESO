####################
Solvers and examples
####################

   .. Common
   
   .. Building solver execs (cmake, spack options)
   .. To build the solver executable, follow the instructions for building NESO in the [top-level README](../../README.md). -->
   
   .. Mesh (re-gen)

   .. In order to generate Nektar++ xml meshes, you'll need `gmsh` and `NekMesh`.
   .. If NESO was installed with spack, then `NekMesh` should already be built.  It can be added to your path with:

   ..    export PATH=$PATH:$(spack location -i nektar%[compiler])/bin

   .. where [compiler] is either 'gcc' or 'oneapi' (both should work if 'spack install' completed without errors.)

   .. Running e.g.s

   .. This script expects to find mpirun on the path and executes with four MPI ranks by default. It looks for a solver executable in the most recently modified spack-build* directory, but this can be overridden using the '-b' option.


.. include:: solver_readmes.md
   :parser: myst_parser.sphinx_


#!/bin/bash
if [ $# -ne 1 ]; then
    echo "usage: cmake_opts [relative_buildpath]"
    exit 1
fi
build_dir=$(realpath "$1")

# Target locations
declare -A target_subdirs
target_subdirs["NESO"]="CMakeFiles/nesolib.dir"
target_subdirs["SimpleSOL_ObjLib"]="solvers/SimpleSOL/CMakeFiles/SimpleSOL_ObjLib.dir"
target_subdirs["SimpleSOL_Exec"]="solvers/SimpleSOL/CMakeFiles/SimpleSOL.dir"
target_subdirs["H3LAPD_ObjLib"]="solvers/H3LAPD/CMakeFiles/H3LAPD_ObjLib.dir"
target_subdirs["H3LAPD_Exec"]="solvers/H3LAPD/CMakeFiles/H3LAPD.dir"

# Grep for build flags, defines
for target_name in ${!target_subdirs[@]}; do
    subdir=${target_subdirs[$target_name]}
    echo
    echo "Target: $target_name"
    echo "  cxx_flags: $(grep CXX_FLAGS $build_dir/$subdir/flags.make|cut -d'=' -f2)"
    echo "  defines: $(grep CXX_DEFINES $build_dir/$subdir/flags.make|cut -d'=' -f2)"
done
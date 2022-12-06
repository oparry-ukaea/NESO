
#===================== Common base for gcc/hipsycl and oneAPI/dpcpp versions ======================
# OS/spack versions hard-coded for now
FROM spack/ubuntu-focal:v0.19.0 as base

ENV spack /opt/spack

COPY NESO-Spack /opt/spack-environment/NESO-Spack
#==================================================================================================

#========================================== gcc/hipsycl ==========================================
FROM base as gcc-hipsycl
#==================================================================================================

#========================================== oneAPI/dpcpp ==========================================
FROM base as oneAPI-dpcpp
# Install intel-oneAPI compilers 
RUN spack install -j 6 intel-oneapi-compilers && spack load intel-oneapi-compilers && spack compiler find
#==================================================================================================

#===================== Common base for gcc/hipsycl and oneAPI/dpcpp versions ======================
# OS/spack versions hard-coded for now
FROM spack/ubuntu-focal:v0.19.0 as base

# Env vars
ENV SYCL_DEVICE_FILTER host
# Variables that are valid in spack ymls (e.g. $spack) are lost because the auto-generated Dockerfile shell-expands them.
# The below is a workaround to replace them:
#ENV spack \$spack
ENV spack /opt/spack

COPY NESO-Spack /opt/spack-environment/NESO-Spack
#==================================================================================================

#========================================== gcc/hipsycl ==========================================
FROM base as gcc-hipsycl
#==================================================================================================

#========================================== oneAPI/dpcpp ==========================================
FROM base as oneAPI-dpcpp

# Env vars 
ENV I_MPI_FABRICS shm

# Install intel-oneAPI compilers 
RUN spack install -j 6 intel-oneapi-compilers && spack load intel-oneapi-compilers && spack compiler find
#==================================================================================================
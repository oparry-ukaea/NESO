# Spack base image + intel-oneapi compilers + NESO-Spack dir


#==============================================================================
# OS/spack versions hard-coded for now
FROM spack/ubuntu-focal:v0.19.0
ENV spack /opt/spack
RUN spack install -j2 intel-oneapi-compilers@2022.1.0 && spack load intel-oneapi-compilers && spack compiler find
COPY NESO-Spack /opt/spack-environment/NESO-Spack
#==============================================================================
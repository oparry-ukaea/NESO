#!/bin/env bash

# Fix these for now
toolchain="gcc-hipsycl"
builder_base_version="0.0.1"
OS_tag="ubuntu:20.04"
num_build_jobs="4"

# Check for prerequisites
if ! command -v docker-compose &> /dev/null
then
    echo "$0: Required command docker-compose not found"
    exit
fi
if ! command -v spack &> /dev/null
then
    echo "$0: Required command spack not found"
    exit
fi

# Set some paths relative to this script
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_root_dir=$(dirname $script_dir)
docker_dir="$repo_root_dir/docker"
#compose_config="$docker_dir/docker-compose.yaml"
spack_root_config="$repo_root_dir/spack.yaml"
spack_config="$docker_dir/spack.yaml"

# Copy root-level config and add container options to it
# Could use an 'include' in the yaml file, but that's not as flexible...
\cp "$spack_root_config" "$spack_config"
printf "  container:\n    images:\n      build: neso/env-builder-base-${toolchain}:${builder_base_version}\n      final: ${OS_tag}\n" >> "$spack_config"
printf "  config:\n    build_jobs: $num_build_jobs" >> "$spack_config"
# Switch to docker dir
cd "$docker_dir"

# Generate Dockerfile
spack containerize > Dockerfile
# Pre-generated Dockerfile doesn't allow for /opt/spack-environment already existing
sed -i 's/mkdir/mkdir\ -p/' Dockerfile
# Use a build cache for spack install
sed -i 's|RUN\ cd /opt/spack-environment\ \&\&\ spack\ env\ activate\ .\ \&\&\ spack\ install|RUN\ --mount=type=cache,target=/opt/spack/var/spack/cache\ cd\ /opt/spack-environment\ \&\&\ spack\ env\ activate\ .\ \&\&\ spack\ install|' Dockerfile


# Build with docker-compose
docker-compose build

# Leave docker dir
cd -
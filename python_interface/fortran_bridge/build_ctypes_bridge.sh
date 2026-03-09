#!/usr/bin/env bash
set -euo pipefail

# Build shared library for the DUDI Python ctypes bridge.
# Does not affect the main Makefile or Fortran-only workflows.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC_DIR="${ROOT_DIR}/src"
BRIDGE_DIR="${ROOT_DIR}/python_interface/fortran_bridge"
PKG_DIR="${ROOT_DIR}/python_interface/dudi"
BUILD_DIR="${BRIDGE_DIR}/build_ctypes"

LIB_OUT="${PKG_DIR}/libpy_dudi_bridge.so"

mkdir -p "${BUILD_DIR}" "${PKG_DIR}"

export LC_ALL=${LC_ALL:-C.UTF-8}
export LANG=${LANG:-C.UTF-8}

FC=${FC:-gfortran}
debug=${DEBUG:-0}

if [[ "$debug" == "1" ]]; then
  echo "• DEBUG build"
  FFLAGS="-O0 -g -fPIC -fopenmp -fcheck=all -finit-real=snan -finit-local-zero -fbacktrace"
  LDFLAGS="-shared -fopenmp -g"
else
  FFLAGS=${FFLAGS:-"-O2 -fPIC -fopenmp"}
  LDFLAGS=${LDFLAGS:-"-shared -fopenmp"}
fi

# DUDI core + batching (order must satisfy .mod dependencies)
FILES_IN_ORDER=(
  "const.f90"
  "comparison_utils.f90"
  "define_types.f90"
  "help.f90"
  "distributions_fun.f90"
  "gu.f90"
  "twobody_fun.f90"
  "inputdata.f90"
  "dataoutmod.f90"
  "integrator.f90"
  "batching.f90"
)

echo "• Building DUDI Fortran core + bridge → ${LIB_OUT}"
MODFLAGS="-J${BUILD_DIR} -I${BUILD_DIR}"

pushd "${SRC_DIR}" >/dev/null
for src in "${FILES_IN_ORDER[@]}"; do
  [[ -f "${src}" ]] || { echo "  (skip) ${src}"; continue; }
  echo "  compile ${src}"
  ${FC} -c ${FFLAGS} ${MODFLAGS} -o "${BUILD_DIR}/${src%.f90}.o" "${src}"
done
popd >/dev/null

echo "  compile py_bridge.f90"
${FC} -c ${FFLAGS} ${MODFLAGS} -o "${BUILD_DIR}/py_bridge.o" "${BRIDGE_DIR}/py_bridge.f90"

echo "  link ${LIB_OUT}"
${FC} ${LDFLAGS} -o "${LIB_OUT}" "${BUILD_DIR}"/*.o

echo "• Built: ${LIB_OUT}"

if ! nm -D "${LIB_OUT}" 2>/dev/null | grep -q " py_dudi_density$"; then
  echo "ERROR: py_dudi_density not found in ${LIB_OUT}"
  exit 1
fi

echo "✓ Bridge build complete."

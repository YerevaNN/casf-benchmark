#!/bin/bash
set -euo pipefail

# Reproduce the isolated Pharmit CLI used by the PLINDER-22 pilot.  The small
# patches include compatibility fixes for current GCC/CMake and the conda-forge
# LEMON package plus an opt-in exact-search mode. The latter disables Pharmit's
# non-conservative higher-order triplet prefilters while retaining independent
# triplet distance bounds and final rigid correspondence validation.

PREFIX="${1:-/mnt/weka/mbedrosian/environments/casf-pharmit-pilot}"
SOURCE_ROOT="${2:-/mnt/weka/mbedrosian/tool-sources/casf-pharmit-pilot}"
PACKAGE_CACHE="${3:-/mnt/weka/mbedrosian/.cache/casf-pharmit-conda}"
MAMBA="${MAMBA:-/opt/conda/bin/mamba}"

mkdir -p "${SOURCE_ROOT}" "${PACKAGE_CACHE}"
CONDA_PKGS_DIRS="${PACKAGE_CACHE}" "${MAMBA}" create -y -p "${PREFIX}" -c conda-forge \
  boost=1.74 cmake=3.27 make jsoncpp eigen lemon gperftools ann curl zlib \
  pkg-config autoconf automake libtool git coincbc

clone_at() {
  local url="$1" destination="$2" revision="$3"
  if [[ ! -d "${destination}/.git" ]]; then
    git clone "${url}" "${destination}"
  fi
  git -C "${destination}" fetch --all --tags
  git -C "${destination}" checkout --detach "${revision}"
}

clone_at https://github.com/FastCGI-Archives/fcgi2.git \
  "${SOURCE_ROOT}/fcgi2" 47f2c03b7771f0ef61d887734ef91e6fa747f837
clone_at https://github.com/openbabel/openbabel.git \
  "${SOURCE_ROOT}/openbabel" cbd4db43f8908b874864280fdc03bf92569eebc1
clone_at https://git.code.sf.net/p/smina/code \
  "${SOURCE_ROOT}/smina-sourceforge" 192d6fdeac509cdee182f978313383390b527961
clone_at https://github.com/tlk00/BitMagic.git \
  "${SOURCE_ROOT}/BitMagic" af3d609ccd873d3c4cb218b14dc4ae5891c38911
clone_at https://github.com/dkoes/pharmit.git \
  "${SOURCE_ROOT}/pharmit" 49d2c0a573f6e58bca62f674b33874bdd42c2274

EXACT_SEARCH_PATCH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pharmit_exact_search.patch"
if git -C "${SOURCE_ROOT}/pharmit" apply --reverse --check "${EXACT_SEARCH_PATCH}" 2>/dev/null; then
  : # already applied
else
  git -C "${SOURCE_ROOT}/pharmit" apply "${EXACT_SEARCH_PATCH}"
fi

pushd "${SOURCE_ROOT}/fcgi2" >/dev/null
autoreconf -i
CC=/usr/bin/gcc CXX=/usr/bin/g++ ./configure --prefix="${PREFIX}"
make -j16
make install
popd >/dev/null

CGICC_TARBALL="${SOURCE_ROOT}/cgicc-3.2.20.tar.gz"
if [[ ! -f "${CGICC_TARBALL}" ]]; then
  "${PREFIX}/bin/curl" -L -o "${CGICC_TARBALL}" \
    https://ftp.gnu.org/gnu/cgicc/cgicc-3.2.20.tar.gz
fi
echo "0f91c5ce8adf406ad35a46152ffb2b9c8d6079604b796d81f48b66d1f068ad4e  ${CGICC_TARBALL}" | sha256sum -c -
if [[ ! -d "${SOURCE_ROOT}/cgicc-3.2.20" ]]; then
  tar -C "${SOURCE_ROOT}" -xzf "${CGICC_TARBALL}"
fi
pushd "${SOURCE_ROOT}/cgicc-3.2.20" >/dev/null
CC=/usr/bin/gcc CXX=/usr/bin/g++ ./configure --prefix="${PREFIX}"
make -j16
make install
popd >/dev/null

if ! grep -q '^#include <ctime>$' "${SOURCE_ROOT}/openbabel/include/openbabel/obutil.h"; then
  sed -i '/#include <math.h>/a #include <ctime>' \
    "${SOURCE_ROOT}/openbabel/include/openbabel/obutil.h"
fi
"${PREFIX}/bin/cmake" --fresh -S "${SOURCE_ROOT}/openbabel" \
  -B "${SOURCE_ROOT}/openbabel/build-casf" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" -DZLIB_LIBRARY="${PREFIX}/lib/libz.so" \
  -DZLIB_INCLUDE_DIR="${PREFIX}/include" -DPYTHON_BINDINGS=OFF -DRUN_SWIG=OFF \
  -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DCMAKE_INSTALL_RPATH="${PREFIX}/lib"
"${PREFIX}/bin/cmake" --build "${SOURCE_ROOT}/openbabel/build-casf" -j16
"${PREFIX}/bin/cmake" --install "${SOURCE_ROOT}/openbabel/build-casf"

"${PREFIX}/bin/cmake" --fresh -S "${SOURCE_ROOT}/smina-sourceforge" \
  -B "${SOURCE_ROOT}/smina-sourceforge/build-casf" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_CXX_STANDARD=14 \
  -DCMAKE_PREFIX_PATH="${PREFIX}" -DBOOST_ROOT="${PREFIX}" \
  -DBoost_NO_SYSTEM_PATHS=ON \
  -DOPENBABEL3_INCLUDE_DIR="${PREFIX}/include/openbabel3" \
  -DOPENBABEL3_LIBRARIES="${PREFIX}/lib/libopenbabel.so" \
  -DEigen3_DIR="${PREFIX}/share/eigen3/cmake" \
  -DCMAKE_BUILD_RPATH="${PREFIX}/lib"
"${PREFIX}/bin/cmake" --build "${SOURCE_ROOT}/smina-sourceforge/build-casf" \
  --target libsmina -j16

# conda-forge LEMON enables several solvers but its umbrella header includes only
# the first one.  Pharmit's packer selects CLP/CBC, so include every enabled API.
python - "${PREFIX}/include/lemon/lp.h" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
s = p.read_text()
s = s.replace('#elif LEMON_HAVE_CPLEX', '#endif\n#ifdef LEMON_HAVE_CPLEX')
s = s.replace('#elif LEMON_HAVE_SOPLEX', '#endif\n#ifdef LEMON_HAVE_SOPLEX')
s = s.replace('#elif LEMON_HAVE_CLP', '#endif\n#ifdef LEMON_HAVE_CLP')
s = s.replace('#elif LEMON_HAVE_CBC', '#endif\n#ifdef LEMON_HAVE_CBC')
p.write_text(s)
PY

rm -rf "${SOURCE_ROOT}/pharmit/src/bm"
cp -a "${SOURCE_ROOT}/BitMagic/src" "${SOURCE_ROOT}/pharmit/src/bm"
if ! grep -q '^#include <chrono>$' "${SOURCE_ROOT}/pharmit/src/SpinLock.h"; then
  sed -i '/#include <cassert>/a #include <chrono>' \
    "${SOURCE_ROOT}/pharmit/src/SpinLock.h"
fi
sed -i 's/target_link_libraries(pharmit jsoncpp /target_link_libraries(pharmit jsoncpp_lib /' \
  "${SOURCE_ROOT}/pharmit/src/CMakeLists.txt"

"${PREFIX}/bin/cmake" --fresh -S "${SOURCE_ROOT}/pharmit/src" \
  -B "${SOURCE_ROOT}/pharmit/src/build-casf" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_CXX_STANDARD=14 \
  -DCMAKE_PREFIX_PATH="${PREFIX}" -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_INSTALL_RPATH="${PREFIX}/lib" -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DBOOST_ROOT="${PREFIX}" -DBoost_NO_SYSTEM_PATHS=ON \
  -DSMINA_DIR="${SOURCE_ROOT}/smina-sourceforge" \
  -DSMINALIB="${SOURCE_ROOT}/smina-sourceforge/build-casf/libsmina.a" \
  -Dbm_SOURCE_DIR="${SOURCE_ROOT}/BitMagic/src" \
  -DLEMON_DIR="${PREFIX}/share/lemon/cmake" \
  -DOpenBabel3_DIR="${PREFIX}/lib/cmake/openbabel3"
"${PREFIX}/bin/cmake" --build "${SOURCE_ROOT}/pharmit/src/build-casf" \
  --target pharmit -j16
"${PREFIX}/bin/cmake" --install "${SOURCE_ROOT}/pharmit/src/build-casf"
"${PREFIX}/bin/pharmit" --help | head -1

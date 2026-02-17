#!/usr/bin/env bash
# =============================================================================
# 03_install_spades.sh
# Workshop: WGS De Novo Assembly — Dr. Omics Lab
# Purpose:  Download and compile SPAdes v4.2.0 from source
#
# Why from source?
#   SPAdes v4.2.0 is not in Ubuntu/Debian apt repositories.
#   We download the tarball from GitHub and compile it ourselves.
#
# Usage: bash scripts/03_install_spades.sh
# Note:  Run from the project root directory.
# =============================================================================

set -euo pipefail

SPADES_VERSION="4.2.0"
TOOLS_DIR="tools"
mkdir -p "${TOOLS_DIR}" logs

echo "======================================================"
echo " Step 3: Installing SPAdes v${SPADES_VERSION} from source"
echo "======================================================"

# Install build dependencies
echo "[1/4] Installing build dependencies..."
sudo apt-get update -qq
sudo apt-get install -y \
    g++ cmake \
    zlib1g-dev \
    build-essential \
    libbz2-dev libz-dev \
    libcurl4-openssl-dev libssl-dev

echo ""
echo "[2/4] Downloading SPAdes v${SPADES_VERSION}..."
cd "${TOOLS_DIR}"
wget -c "https://github.com/ablab/spades/archive/refs/tags/v${SPADES_VERSION}.tar.gz"

echo ""
echo "[3/4] Extracting and compiling..."
tar -xzf "v${SPADES_VERSION}.tar.gz"
cd "spades-${SPADES_VERSION}"
./spades_compile.sh 2>&1 | tee "../../logs/03_spades_compile.log"

echo ""
echo "[4/4] Verifying installation..."
SPADES_BIN="$(pwd)/bin/spades.py"
if "${SPADES_BIN}" --version &>/dev/null; then
    echo "✓ SPAdes installed successfully!"
    echo "  Executable: ${SPADES_BIN}"
else
    echo "ERROR: Compilation may have failed. Check logs/03_spades_compile.log"
    exit 1
fi

echo ""
echo "======================================================"
echo " Installation complete."
echo " To run SPAdes, use the full path:"
echo "   $(realpath "${SPADES_BIN}")"
echo ""
echo " Or open WSL directly from:"
echo "   ${TOOLS_DIR}/spades-${SPADES_VERSION}/bin/"
echo " and run: ./spades.py"
echo "======================================================"

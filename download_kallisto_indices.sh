#!/bin/bash

# Script to download Kallisto indices from the official repository
# Usage: ./download_kallisto_indices.sh [species] [type]
# species: all, human, mouse, dog, monkey, zebrafish (default: all)
# type: standard, nac (default: standard)

set -e

SPECIES="${1:-all}"
TYPE="${2:-standard}"
BASE_URL="https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1"
DATA_DIR="data/kallisto_indices"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to download and extract index
download_index() {
    local species=$1
    local type=$2
    local filename="${species}_index_${type}.tar.xz"
    local url="${BASE_URL}/${filename}"
    local target_dir="${DATA_DIR}/${species}"
    
    echo -e "${YELLOW}Downloading ${species} ${type} index...${NC}"
    
    # Create target directory
    mkdir -p "${target_dir}"
    
    # Download the file
    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "/tmp/${filename}" "${url}"
    elif command -v curl &> /dev/null; then
        curl -L -# -o "/tmp/${filename}" "${url}"
    else
        echo -e "${RED}Error: Neither wget nor curl is installed${NC}"
        exit 1
    fi
    
    # Extract the archive
    echo -e "${YELLOW}Extracting ${species} index...${NC}"
    tar -xf "/tmp/${filename}" -C "${target_dir}"
    
    # Clean up
    rm "/tmp/${filename}"
    
    echo -e "${GREEN}✓ ${species} ${type} index downloaded successfully${NC}"
}

# Function to check if index already exists
check_existing() {
    local species=$1
    local target_dir="${DATA_DIR}/${species}"
    
    if [ -f "${target_dir}/index.idx" ] && [ -f "${target_dir}/t2g.txt" ]; then
        return 0
    else
        return 1
    fi
}

# Main logic
echo -e "${GREEN}Kallisto Index Downloader${NC}"
echo "========================="

# Create base data directory if it doesn't exist
mkdir -p "${DATA_DIR}"

# Available species
declare -a species_list=("human" "mouse" "dog" "monkey" "zebrafish")

# Download based on user selection
if [ "${SPECIES}" = "all" ]; then
    echo -e "${YELLOW}Downloading all species indices (type: ${TYPE})${NC}"
    for sp in "${species_list[@]}"; do
        if check_existing "${sp}"; then
            echo -e "${GREEN}✓ ${sp} index already exists, skipping...${NC}"
        else
            download_index "${sp}" "${TYPE}"
        fi
    done
elif [[ " ${species_list[@]} " =~ " ${SPECIES} " ]]; then
    if check_existing "${SPECIES}"; then
        echo -e "${YELLOW}Index for ${SPECIES} already exists.${NC}"
        read -p "Do you want to re-download? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            download_index "${SPECIES}" "${TYPE}"
        else
            echo -e "${GREEN}Skipping download.${NC}"
        fi
    else
        download_index "${SPECIES}" "${TYPE}"
    fi
else
    echo -e "${RED}Error: Invalid species '${SPECIES}'${NC}"
    echo "Available species: ${species_list[*]}"
    echo "Usage: $0 [species] [type]"
    echo "  species: all, human, mouse, dog, monkey, zebrafish (default: all)"
    echo "  type: standard, nac (default: standard)"
    exit 1
fi

# Verify downloaded files
echo -e "\n${GREEN}Verification:${NC}"
for sp in "${species_list[@]}"; do
    if [ -d "${DATA_DIR}/${sp}" ]; then
        if [ -f "${DATA_DIR}/${sp}/index.idx" ] && [ -f "${DATA_DIR}/${sp}/t2g.txt" ]; then
            echo -e "  ${GREEN}✓${NC} ${sp}: index.idx and t2g.txt present"
        else
            echo -e "  ${RED}✗${NC} ${sp}: missing files"
        fi
    fi
done

echo -e "\n${GREEN}Download complete!${NC}"
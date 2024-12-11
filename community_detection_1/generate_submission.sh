#!/bin/bash

# generate directories if they do not exist
if [ ! -d "./networks" ]; then mkdir "./networks"; fi
if [ ! -d "./results/subchallenge1" ]; then mkdir -p "./results/subchallenge1"; fi
if [ ! -d "./results/subchallenge2" ]; then mkdir -p "./results/subchallenge2"; fi


# Download networks from Synapse
./scripts/fetch_networks.py

# Subchallenge 1
# network 1 - protein-protein interaction
./scripts/sc1_n1_ppi.py

# network 2 - protein-protein interaction
./scripts/sc1_n2_ppi.sh

# network 3 - signal directed
./scripts/sc1_n3_signal_directed.sh

# network 4 - coexpression
./scripts/sc1_n4_coexpr.py

# network 5 - cancer
./scripts/sc1_n5_cancer.sh

# network 5 - homology
./scripts/sc1_n6_homology.sh

# Subchallenge 2
./scripts/sc2_overlap.py

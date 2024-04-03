#!/usr/bin/env bash

# ~ Load the virtual environment
SD=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# ==============================================================================
# MERRIN INPUTS
# ==============================================================================
SBML="${SD}/ecoli-small/metabolic_network.sbml" # Metabolic network in SBML
OBJ="Growth"                                    # Objective reaction
PKN="${SD}/ecoli-small/pkn.txt"                 # Prior Knwoledge Network file
# Observations description file in JSON format:
# `timeseries_kft.json`: timeseries observation with kinetics, fluxomics and
#                        transcriptomics data
# `timeseries_t.json`: timeseries observation with transcriptomics data only
OBSERVATIONS="${SD}/ecoli-small/timeseries_kft.json"

# ==============================================================================
# MERRIN INPUTS -- OPTIONAL
# ==============================================================================
OPTIMISATION="subsetmin"      # Optimisation mode, either `all` or `subsetmin`
PROJECTION="network"    # Projection mode, either `network` or `node`
OUTPUT_CSV="${SD}/merrin-ecoli-small.csv"

# ==============================================================================
# Start MERRIN
# ==============================================================================
merrin \
    -sbml ${SBML} -pkn ${PKN} -obj ${OBJ} -obs ${OBSERVATIONS} \
    --projection ${PROJECTION} --optimisation ${OPTIMISATION} \
    -out ${OUTPUT_CSV}
echo "Regulatory rules have been saved in: ${OUTPUT_CSV}"

#!/usr/bin/env bash

# ~ Load the virtual environment
SD=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# ==============================================================================
# INSTANCE
# ==============================================================================
# Instance among:
# - toy, with the objective reaction `Growth`
# - core-regulated, with the OBJective function `Growth`
# - large-scale, with the OBJective function `VGRO`
INSTANCE="core-regulated"   # Instance
OBJ="Growth"                # Objective reaction

# ==============================================================================
# MERRIN INPUTS
# ==============================================================================
SBML="${SD}/instances/${INSTANCE}/metabolic_network.sbml"
PKN="${SD}/instances/${INSTANCE}/pkn.txt"
OBSERVATIONS="${SD}/instances/${INSTANCE}/timeseries_kft.json"

# ==============================================================================
# MERRIN INPUTS -- OPTIONAL
# ==============================================================================
OPTIMISATION="subsetmin" # Optimisation mode, either `all` or `subsetmin`
PROJECTION="trace"       # Projection mode, either `network`, `node`, or `trace`
OUTPUT_CSV="${SD}/merrin-${INSTANCE}.${OPTIMISATION}.${PROJECTION}.csv"

# ==============================================================================
# Start MERRIN
# ==============================================================================
echo "Regulatory rules will been saved in: ${OUTPUT_CSV}"
merrin \
    -sbml ${SBML} -pkn ${PKN} -obj ${OBJ} -obs ${OBSERVATIONS} \
    --projection ${PROJECTION} --optimisation ${OPTIMISATION} \
        | tee ${OUTPUT_CSV}

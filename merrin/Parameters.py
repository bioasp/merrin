#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

from enum import Enum

#ENUM###########################################################################

class UpdateMode(Enum):
    SYNCHRONOUS = 'synchronous'
    ASYNCHRONOUS = 'asynchronous'

#DATA#PROCESSING################################################################

THRESHOLD = 10**-6
SIMULATION_TIMESTEP = 0.01
BIOMASS = 'biomass'

#DATA#SOURCE####################################################################

DATA_TYPE_COLUMN = 'tag'
PROTEOMIC_DATA = 'proteomic'
FLUXOMIC_DATA = 'fluxomic'
CINETIC_DATA = 'cinetic'

#DATA#ANNOTATION################################################################

ANNOTATION_COLUMN = 'annotation'
SUBSTRATE_DEPLETION = 'substrate depletion'
REGULATORY_UPDATE   = 'regulatory update'

#RESOLUTION#####################################################################

MAX_OBS_BUFFER = 0

DATA_ERROR_EPSILON = 0.1
LP_MODEL_EPSILON = 10**-4

UPDATE_MODE = UpdateMode.SYNCHRONOUS

#RESULTS#OUTPUT#################################################################

OUTPUT_DIR = './out'

# ==============================================================================
# Imports
# ==============================================================================
from os import path

from merrin.asp.instantiater import (
    instantiate_parameters,
    instantiate_networks,
    instantiate_mn,
    instantiate_pkn,
    instantiate_observations,
    instantiate_trace
)

# ==============================================================================
# Globals
# ==============================================================================
SD: str = path.dirname(path.abspath(__file__))
ASP_MODEL_LEARN: str = f'{SD}/model/learn.lp'
ASP_MODEL_LEARN_FROM_TRACE: str = f'{SD}/model/learn_from_trace.lp'

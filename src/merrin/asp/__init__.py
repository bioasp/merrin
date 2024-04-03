# ==============================================================================
# Imports
# ==============================================================================
from os import path

from merrin.asp.instantiater import (
    instantiate_parameters,
    instantiate_networks,
    instantiate_observations
)

# ==============================================================================
# Globals
# ==============================================================================
SD: str = path.dirname(path.abspath(__file__))
ASP_MODEL_LEARN: str = f'{SD}/model/learn.lp'
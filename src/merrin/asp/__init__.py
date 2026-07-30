# ==============================================================================
# Imports
# ==============================================================================
from os import path

from merrin.asp.instantiater import (
    instantiate_mn,  # noqa
    instantiate_networks,  # noqa
    instantiate_observations,  # noqa
    instantiate_parameters,  # noqa
    instantiate_pkn,  # noqa
    instantiate_trace_domain,  # noqa
)

# ==============================================================================
# Globals
# ==============================================================================
SD: str = path.dirname(path.abspath(__file__))
ASP_MODEL_LEARN: str = f'{SD}/model/learn.lp'
ASP_MODEL_LEARN_FROM_TRACE: str = f'{SD}/model/learn_from_trace.lp'

# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from typing import Literal

from os.path import dirname
from json import load
from enum import Enum
from pandas import DataFrame, read_csv


# ==============================================================================
# Type Association
# ==============================================================================
Datatype = Literal['fluxomics', 'kinetics', 'transcriptomics']


# ==============================================================================
# Timeseries object
# ==============================================================================
class Observation:

    class Datatype(Enum):
        FLUXOMICS = 'fluxomics'
        KINETICS = 'kinetics'
        TRANSCRIPTOMICS = 'transcriptomics'

    def __init__(self: Observation) -> None:
        self.id: str | None = None
        self.types: set[Datatype] = set()
        self.cstr_mutations: dict[str, bool] = {}
        self.cstr_bounds: dict[str, tuple[float, float]] = {}
        self.data: DataFrame | None = None

    # ==========================================================================
    # Data clean up
    # ==========================================================================
    @classmethod
    def __data_binarized(cls: type[Observation], data: DataFrame) -> DataFrame:
        raise NotImplementedError()

    @classmethod
    def __data_remove_redundancy(cls: type[Observation], data: DataFrame) \
            -> DataFrame:
        raise NotImplementedError()

    # ==========================================================================
    # Load from json files
    # ==========================================================================
    @classmethod
    def load_json(cls: type[Observation], ts_json: str) -> list[Observation]:
        observations: list[Observation] = []
        with open(ts_json, 'r', encoding='utf-8') as file:
            instance_dir: str = dirname(ts_json)
            instance_description: dict = load(file)
            for i, ts in enumerate(instance_description):
                # --------------------------------------------------------------
                # Retrieve values
                # --------------------------------------------------------------
                # FIXME add a checkup on the input datatypes
                ident: str = str(i)
                types: set[Datatype] = {dt.lower() for dt in ts['type']}
                mutations: dict[str, bool] = {
                    n: value == 1
                    for n, value in ts['constraints']['mutations'].items()
                }
                bounds: dict[str, tuple[float, float]] = {
                    r: (float(lb), float(ub))
                    for r, (lb, ub) in ts['constraints']['bounds'].items()
                }
                # --------------------------------------------------------------
                # Load and clean up the DataFrame
                # --------------------------------------------------------------
                csv_file: str = instance_dir + '/' + ts['file']
                observations_df: DataFrame = \
                    read_csv(csv_file).set_index('Time')
                # --------------------------------------------------------------
                # Retrieve values
                # --------------------------------------------------------------
                observation: Observation = Observation()
                observation.id = ident
                observation.types = types
                observation.cstr_mutations = mutations
                observation.cstr_bounds = bounds
                observation.data = observations_df
                observations.append(observation)
        return observations

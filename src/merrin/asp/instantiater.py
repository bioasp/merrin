# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from typing import Iterable
from numpy import isnan
from pandas import DataFrame, Series

import merrin.asp.template as Template
from merrin.datastructure import MetabolicNetwork, Observation


# ==============================================================================
# Metabolic and Prior Knowledge Networks
# ==============================================================================
def instantiate_networks(mn: MetabolicNetwork, objective: str,
                         pkn: Iterable[tuple[str, int, str]],
                         max_clause: int = 20) -> list[str]:
    # --------------------------------------------------------------------------
    # Instantiate ASP constraints for the metabolic network
    # --------------------------------------------------------------------------
    assert mn.irreversible
    mn_asp: list[str] = __instantiate_mn(mn, objective)
    # --------------------------------------------------------------------------
    # Instantiate ASP constraints for the PKN
    # --------------------------------------------------------------------------
    # ~ Instantiate ASP constraints
    pkn_asp: list[str] = __instantiate_pkn(pkn, max_clause=max_clause)
    # --------------------------------------------------------------------------
    # Return
    # --------------------------------------------------------------------------
    return mn_asp + pkn_asp


def __instantiate_mn(mn: MetabolicNetwork, objective: str) -> list[str]:
    mn_asp: list[str] = [Template.header('Metabolic Network', '=')]
    # --------------------------------------------------------------------------
    # Objective
    # --------------------------------------------------------------------------
    mn_asp.append(Template.header('Objective', '-'))
    mn_asp.append(Template.MetabolicNetwork.objective(objective))
    # --------------------------------------------------------------------------
    # External Metabolites
    # --------------------------------------------------------------------------
    mn_asp.append(Template.header('External metabolite', '-'))
    mn_asp.extend(
        Template.MetabolicNetwork.metabolite_external(m_ext)
        for m_ext in sorted(mn.metabolites(external=True, internal=False))
    )
    # --------------------------------------------------------------------------
    # Reactions
    # --------------------------------------------------------------------------
    mn_asp.append(Template.header('Reactions', '-'))
    stoichiometry: dict[tuple[str, str], float] = mn.stoichiometry()
    for r in sorted(mn.reactions()):
        reactant_asp: list[str] = [
            Template.MetabolicNetwork.reactant(r, m, stoichiometry[(m, r)])
            for m, r_ in stoichiometry
            if r_ == r and stoichiometry[(m, r)] < 0
        ]
        product_asp: list[str] = [
            Template.MetabolicNetwork.product(
                r, m, stoichiometry[(m, r)])
            for m, r_ in stoichiometry
            if r_ == r and stoichiometry[(m, r)] > 0
        ]
        lb, ub = mn.bound(r)
        assert 0 <= lb <= ub
        mn_asp.append(Template.comment(r))
        mn_asp.append(Template.MetabolicNetwork.bounds(r, lb, ub))
        if lb == ub:
            mn_asp.append(Template.MetabolicNetwork.fixed_flux(r))
        mn_asp.extend(sorted(reactant_asp))
        mn_asp.extend(sorted(product_asp))
    # --------------------------------------------------------------------------
    # Reversible reactions
    # --------------------------------------------------------------------------
    reversible_reactions: list[tuple[str, str]] = \
        list(mn.previously_reversible_reactions().values())
    if len(reversible_reactions) != 0:
        mn_asp.append(Template.header('Reversible Reactions', '-'))
        for (rf, rr) in reversible_reactions:
            mn_asp.append(Template.MetabolicNetwork.reversible(rf, rr))
    # --------------------------------------------------------------------------
    # Return
    # --------------------------------------------------------------------------
    return mn_asp


def __instantiate_pkn(pkn: Iterable[tuple[str, int, str]],
                      max_clause: int = 20) -> list[str]:
    pkn_asp: list[str] = [Template.header('Prior Knowledge Network', '=')]
    # Convert the PKN to a more convenient form
    formatted_pkn: dict[str, list[tuple[str, int]]] = {}
    for u, s, v in pkn:
        formatted_pkn.setdefault(v, []).append((u, s))
    # --------------------------------------------------------------------------
    # Nodes
    # --------------------------------------------------------------------------
    for n, in_edges in formatted_pkn.items():
        pkn_asp.append(Template.comment(n))
        pkn_asp.append(Template.PKN.node(n, fixed=False))
        pkn_asp.append(
            Template.PKN.max_clause(n, len(in_edges), max_clause=max_clause)
        )
        for u, s in sorted(in_edges):
            pkn_asp.append(Template.PKN.edge(u, n, s))

    # --------------------------------------------------------------------------
    # Return
    # --------------------------------------------------------------------------
    return pkn_asp


# ==============================================================================
# Observations
# ==============================================================================
def instantiate_observations(mn: MetabolicNetwork, objective: str,
                             observations: Iterable[Observation],
                             epsilon: float = 10**-9) -> list[str]:
    # --------------------------------------------------------------------------
    # Instantiate ASP constraints for the PKN
    # --------------------------------------------------------------------------
    observations_asp: list[str] = []
    for observation in observations:
        observations_asp.extend(
            __instantiate_observation(
                mn, observation, objective, epsilon=epsilon
            )
        )
    # --------------------------------------------------------------------------
    # Return
    # --------------------------------------------------------------------------
    return observations_asp


def __instantiate_observation(mn: MetabolicNetwork, observation: Observation,
                              objective: str, epsilon: float = 10**-9) \
        -> list[str]:
    observation_id: str = observation.id  # type: ignore
    observation_asp: list[str] = [
        Template.header(f'Observation <{observation_id}>', '=')
    ]
    # --------------------------------------------------------------------------
    # Special constraints: mutation
    # --------------------------------------------------------------------------
    if len(observation.cstr_mutations) != 0:
        observation_asp.append(
            Template.header('Mutations', '-')
        )
        observation_asp.extend(
            Template.TimeSeries.constraint_mutant(
                observation_id, n_, v
            )
            for n, v in observation.cstr_mutations.items()
            for n_ in mn.previously_reversible_reactions().get(n, [n])
        )
    # --------------------------------------------------------------------------
    # Convert observation data according to the irreversible metabolic network
    # --------------------------------------------------------------------------
    data: DataFrame | None = observation.data
    assert data is not None
    data = data.copy()
    for r, (rf, rr) in mn.previously_reversible_reactions().items():
        if r in data.columns:
            data[rf] = data.loc[:, r].apply(lambda x:  x if x > 0 else 0)
            data[rr] = data.loc[:, r].apply(lambda x: -x if x < 0 else 0)
            data.drop(columns=[r])
    # --------------------------------------------------------------------------
    # Observation of each timestep
    # --------------------------------------------------------------------------
    previous_time: int | None = None
    for time in data.index:
        timestep: Series = data.loc[time]
        if previous_time is not None:
            observation_asp.extend(
                __instantiate_bounds(
                    (observation_id, previous_time), observation,
                    (data.loc[previous_time], timestep), mn, tau=0.01
                )
            )
            observation_asp.append(
                Template.TimeSeries.observation_transition(
                    (observation_id, previous_time), (observation_id, time)
                )
            )
        observation_asp.extend(
            __instantiate_timestep(
                (observation_id, time), timestep,
                observation.types, objective, epsilon  # type: ignore
            )
        )
        previous_time = time
    if previous_time is not None:
        observation_asp.extend(
            __instantiate_bounds(
                (observation_id, previous_time), observation,
                (data.loc[previous_time], data.loc[previous_time]),
                mn, tau=0.01
            )
        )
    return observation_asp


def __instantiate_timestep(ident: tuple[str, int], timestep: Series,
                           datatypes: set[str], objective: str,
                           epsilon: float = 10**-9) -> list[str]:
    timestep_asp: list[str] = [Template.header(f'Time: {ident[1]}', '-')]
    if objective in timestep.keys() and \
            ('fluxomics' in datatypes or 'kinetics' in datatypes):
        timestep_asp.append(Template.comment('Growth'))
        timestep_asp.append(Template.TimeSeries.observation_growth(
            ident, float(timestep[objective])
        ))
    timestep_asp.append(Template.comment('Binarized Observations'))
    for element in timestep.keys():
        if isnan(timestep[element]):
            continue
        timestep_asp.append(
            Template.TimeSeries.observation_value(
                ident, element, timestep[element], epsilon
            )
        )
    return timestep_asp


def __instantiate_bounds(ident: tuple[str, int], observation: Observation,
                         timesteps: tuple[Series, Series],
                         mn: MetabolicNetwork, tau: float = 0.01) -> list[str]:
    datatypes: set[str] = observation.types  # type: ignore
    bounds_asp: list[str] = []
    if 'fluxomics' not in datatypes and 'kinetics' not in datatypes:
        return bounds_asp
    bounds_asp.append(Template.comment('Flux constraints'))
    for m in mn.metabolites(internal=False):
        r: str | None = mn.exchange(m)
        assert r is not None
        if mn.stoichiometry()[(m, r)] > 0:
            continue
        lb, ub = observation.cstr_bounds.get(r, mn.bound(r))
        if 'kinetics' in datatypes:
            assert 'biomass' in timesteps[0].keys()
            w: float = timesteps[0][m]
            biomass: float = timesteps[0]['biomass']
            uptake: float = max(0, __compute_uptake(w, biomass, tau=tau))
            if uptake < ub:
                bounds_asp.append(Template.TimeSeries.observation_exchange(
                    ident, r, lb, uptake
                ))
        elif 'fluxomics' in datatypes and \
                timesteps[0][m] > 0 and timesteps[1][m] == 0:
            if timesteps[0][r] < ub:
                bounds_asp.append(Template.TimeSeries.observation_exchange(
                    ident, r, lb, timesteps[0][r]
                ))
    return bounds_asp


def __compute_uptake(w: float, biomass: float, tau: float = 0.01) -> float:
    return w / (biomass * tau)


# ==============================================================================
# Parameters
# ==============================================================================
def instantiate_parameters(max_gap: int, max_error: float,
                           epsilon: float) -> list[str]:
    parameters_asp: list[str] = [
        Template.header('Parameters', '='),
        Template.Parameter.max_gap(max_gap),
        Template.Parameter.max_error(max_error),
        Template.Parameter.lp_epsilon(epsilon)
    ]
    return parameters_asp

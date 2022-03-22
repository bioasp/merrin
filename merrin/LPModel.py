#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

from pulp import (
    COIN_CMD,
    PULP_CBC_CMD,
    LpStatus,
    LpMaximize,
    LpProblem,
    LpVariable,
    LpBinary,
    lpSum,
)

#AUXILIARY#FUNCTIONS##MODELS####################################################

LP_SOLVER = COIN_CMD if COIN_CMD().available() else PULP_CBC_CMD


def get_base_fba_lp_model(time, internal_metabolites, reactions, bounds, stoichiometry, obj, removeObj=False):
    # ------------------------------------------------------------------------ #
    # Linear problem with minimization
    # ------------------------------------------------------------------------ #
    fba = LpProblem('FBA', LpMaximize)

    # ------------------------------------------------------------------------ #
    # The variables
    # ------------------------------------------------------------------------ #
    f = {(t, r): LpVariable(f'f_({t}_{r})', lowBound=bounds[(t, r)][0], upBound=bounds[(t, r)][1])
         for t in time for r in reactions}

    # ------------------------------------------------------------------------ #
    # The objective function
    # ------------------------------------------------------------------------ #
    # Maximise the objective reaction
    if not removeObj:
        fba += lpSum(f[(t, obj)]
                     for t in time if obj in reactions), 'Objective_reaction'

    # ------------------------------------------------------------------------ #
    # The constraints
    # ------------------------------------------------------------------------ #
    # Steady state: \forall m \in M_internal, \sum_{r \in R} S_{r,m} * f_r = 0
    for t in time:
        for m in internal_metabolites:
            fba += \
                lpSum(f[(t, r)] * stoichiometry.get((r, m), 0) for r in reactions) == 0, \
                f'steady_state_{t}_{m}'

    return fba, f


def add_inhibition_constraints(fba, f, reactions, reactions_states):
    # ------------------------------------------------------------------------ #
    # The constraints
    # ------------------------------------------------------------------------ #
    #
    for ((t1, t2), r, s) in reactions_states:
        if r in reactions and s == -1:
            t = (t1, t2)
            fba += f[(t, r)] == 0, f'inhibition_{t}_{r}'

    return fba


def add_forced_state(fba, f, reactions, reactions_states, epsilon=10**-5):
    # ------------------------------------------------------------------------ #
    # The constraints
    # ------------------------------------------------------------------------ #
    #
    for ((t1, t2), r, s) in reactions_states:
        t = (t1, t2)
        if r in reactions and s == 1:
            fba += f[(t, r)] >= epsilon, f'activated_{t}_{r}'
        if r in reactions and s == -1:
            fba += f[(t, r)] == 0, f'activated_{t}_{r}'

    return fba


def add_maximised_forced_state(fba, f, reactions, reactions_states, bounds, epsilon=10**-5):
    # ------------------------------------------------------------------------ #
    # The parameters
    # ------------------------------------------------------------------------ #
    time = set([t for (t, _, _) in reactions_states])

    # ------------------------------------------------------------------------ #
    # The variables
    # ------------------------------------------------------------------------ #
    y = {t: LpVariable(f'y_({t})', cat=LpBinary) for t in time}

    # ------------------------------------------------------------------------ #
    # The objective function
    # ------------------------------------------------------------------------ #
    #
    fba += lpSum(y[t] for t in time), 'Objective function'

    # ------------------------------------------------------------------------ #
    # The constraints
    # ------------------------------------------------------------------------ #
    #
    for ((t1, t2), r, s) in reactions_states:
        t = (t1, t2)
        if r in reactions and s == 1:
            fba += f[(t, r)] >= epsilon * y[t], f'activated_{t}_{r}'
        if r in reactions and s == -1:
            fba += f[(t, r)] <= (1 - y[t]) * \
                bounds[(t, r)][1], f'activated_{t}_{r}'

    return fba, y


def solve_rfba(fba_lp_model, lp_variables, obj_reaction):
    # Coin Branch and Cut solver is used to solve the instanced model
    solver = LP_SOLVER(msg=0, timeLimit=5)

    fba_lp_model.solve(solver)

    # Retrieve the values of the objective reaction for each time step
    obj_values = {}
    for (t, r) in lp_variables.keys():
        if r == obj_reaction:
            r_value = lp_variables[(t, r)].varValue
            obj_values[t] = r_value

    solver_status = LpStatus[fba_lp_model.status]

    return solver_status, obj_values


def solve_lp_model(lp_model, lp_variables):
    # Coin Branch and Cut solver is used to solve the instanced model
    solver = LP_SOLVER(msg=0)

    lp_model.solve(solver)

    solver_status = LpStatus[lp_model.status]

    # Retrieve the values of the objective reaction for each time step
    values = ([], [])
    for rn in lp_variables.keys():
        y_value = lp_variables[rn].varValue
        if y_value == 1:
            values[0].append(rn)
        else:
            values[1].append(rn)

    return solver_status, values


#CODE###########################################################################

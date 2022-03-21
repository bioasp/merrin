#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

from .LPModel import *
from .Parameters import (
    THRESHOLD,
    DATA_ERROR_EPSILON,
    BIOMASS
)

from .utils import count_and_time


#CLASS#PROPAGATOR###############################################################


#==============================================================================#
# FBA PROPAGATOR
#==============================================================================#

class FBAPropagator():

    #--------------------------------------------------------------------------#
    # Initialisation
    #--------------------------------------------------------------------------#

    def __init__(self, mn, pkn, obs, obj, data_error=DATA_ERROR_EPSILON):
        self._fluxprops = []
        self._mn = mn
        self._pkn = pkn
        self._obs = obs
        self._obj = obj
        self.__data_error = data_error

        # Memory
        self.generated_nogoods = []
        self.to_add_nogoods = []

    def init(self, init):
        self.reset()
        for _ in range((init.number_of_threads)):
            fluxprop = FBAChecker(
                self._mn,
                self._pkn,
                self._obs,
                self._obj,
                data_error=self.__data_error
            )
            fluxprop.init(init)
            self._fluxprops.append(fluxprop)

    #--------------------------------------------------------------------------#
    # Memory management for consecutive resolutions
    #--------------------------------------------------------------------------#

    def reset(self):
        self._fluxprops = []

    def reset_memory(self):
        self.generated_nogoods = []

    #--------------------------------------------------------------------------#
    # Propagator's functions
    #--------------------------------------------------------------------------#

    def propagate(self, control, change):
        fluxprop = self._fluxprops[control.thread_id]

        #......................................................................#
        # Update the variable assignment
        #......................................................................#
        fluxprop.propagate(change)

        #......................................................................#
        # Check if the assignment respects the linear constraints
        #......................................................................#
        status, nogoods = fluxprop.propagate_verification(change)

        #......................................................................#
        # Add a nogood if it exists
        #......................................................................#
        if status and len(nogoods) != 0:
            self.to_add_nogoods.extend(nogoods)
        while len(self.to_add_nogoods) != 0:
            nogood = self.to_add_nogoods.pop()
            self.generated_nogoods.append(fluxprop.nogood2str(nogood))
            # print(fluxprop.nogood2str(nogood))
            if not control.add_nogood(nogood):
                return

    def undo(self, thread_id, _, change):
        fluxprop = self._fluxprops[thread_id]

        #......................................................................#
        # Update the variable assignment
        #......................................................................#
        fluxprop.undo(change)

    def check(self, control):
        fluxprop = self._fluxprops[control.thread_id]

        #......................................................................#
        # Check if the assignment respects the linear constraints
        #......................................................................#
        status, nogoods = fluxprop.check_verification(control)

        #......................................................................#
        # Add a nogood if it exists
        #......................................................................#
        if status and len(nogoods) != 0:
            self.to_add_nogoods.extend(nogoods)
        while len(self.to_add_nogoods) != 0:
            nogood = self.to_add_nogoods.pop()
            self.generated_nogoods.append(fluxprop.nogood2str(nogood))
            # print(fluxprop.nogood2str(nogood))
            if not control.add_nogood(nogood):
                return

#==============================================================================#
# FBA CHECKER
#==============================================================================#


class FBAChecker:

    #--------------------------------------------------------------------------#
    # Initialisation
    #--------------------------------------------------------------------------#

    def __init__(self, mn, pkn, obs, obj, data_error=DATA_ERROR_EPSILON):
        """
        """
        #......................................................................#
        # Variable assignment
        #......................................................................#
        self._map_p2l = {}
        self._map_l2p = {}

        self._mss_states = set()
        self._reg_states = set()
        self._node_states = set()
        self._clause_states = set()

        #......................................................................#
        # Generated nogoods
        #......................................................................#
        self.nogoods = []

        #......................................................................#
        # Error rate
        #......................................................................#
        self.__data_error = data_error

        #......................................................................#
        # Observations
        #......................................................................#
        self.__observations = obs
        self.__time_obs = obs.index.values.tolist()
        self.__time_model = set()

        #......................................................................#
        # Regulatory metabolic network informations
        #......................................................................#
        self.__metabolic_network = mn
        self.__obj = obj
        self.__regulators = set(pkn.nodes)      \
            .difference(mn.get_reactions())     \
            .difference(mn.get_inputs())        \
            .difference(mn.get_outputs())
        self.__node_candidates = set(pkn.nodes) \
            .difference(self.__regulators)      \
            .intersection(mn.get_reactions())

        #......................................................................#
        # Metabolic bounds computations
        #......................................................................#
        self.__bounds = {}

        # For import reactions
        for t in self.__time_obs:
            for r, m in mn.get_input_reactions():
                if BIOMASS in obs.columns and obs.loc[t][BIOMASS] is not None:
                    biomass = obs.loc[t][BIOMASS]
                    C = obs.loc[t][m]
                    r_bound = mn.get_transport_reaction_bounds(
                        r, m, C, biomass)
                    self.__bounds[(t, r)] = r_bound
                elif self.__observations.loc[t][m] is not None:
                    C = self.__observations.loc[t][m]
                    if C == 0:
                        self.__bounds[(t, r)] = (0, 0)
                    else:
                        self.__bounds[(t, r)] = mn.get_reaction_bounds(r)
                else:
                    self.__bounds[(t, r)] = mn.get_reaction_bounds(r)

        # For internal reactions
        internal_reactions = mn.get_reactions() \
            .difference(r for r, _ in mn.get_input_reactions())
        for t in self.__time_obs:
            for r in internal_reactions:
                self.__bounds[(t, r)] = \
                    self.__metabolic_network.get_reaction_bounds(r)

        #......................................................................#
        # Extract objective values
        #......................................................................#
        self.__objective_values = {}
        for t in self.__time_obs:
            if self.__obj in obs.columns:
                self.__objective_values[t] = obs.loc[t][obj]

        #......................................................................#
        # Cache
        #......................................................................#
        self._cache_valid = set()

    def init(self, init) -> None:
        """
        """
        #......................................................................#
        # Initialise the system state
        #......................................................................#
        self._map_p2l = {'clause': {}, 'node': {}, 'w': {}, 'v': {}}
        self._map_l2p = {}

        #......................................................................#
        # Define 'to watch' variables
        #......................................................................#
        reactions = self.__metabolic_network.get_reactions()
        candidates = self.__node_candidates
        inputs = self.__metabolic_network.get_inputs()

        for atom in init.symbolic_atoms:
            guess_lit = init.solver_literal(atom.literal)
            guess_truth = init.assignment.value(guess_lit)

            symbol = atom.symbol

            if guess_truth is False:
                continue

            #..................................................................#
            # 'w_obs(t, r, v)' - regulation state
            #..................................................................#
            if symbol.name == 'w_obs':
                t = (symbol.arguments[0].arguments[0].number,
                     symbol.arguments[0].arguments[1].number)
                r = symbol.arguments[1].string
                v = symbol.arguments[2].number
                if r in reactions and t in self.__time_obs:
                    self._map_p2l['w'].setdefault(t, {})[(r, v)] = guess_lit
                    self._map_l2p[guess_lit] = ('w', t, r, v)
                    init.add_watch(guess_lit)
                    if guess_truth is True:
                        self._reg_states.add((t, r, v))

            #..................................................................#
            # 'clause(n, c, a, v)' - clause
            #..................................................................#
            elif symbol.name == 'clause':
                n = symbol.arguments[0].string
                c = symbol.arguments[1].number
                a = symbol.arguments[2].string
                v = symbol.arguments[3].number
                self._map_p2l['clause'][(n, c, a, v)] = guess_lit
                self._map_l2p[guess_lit] = ('clause', n, c, a, v)
                init.add_watch(guess_lit)
                if guess_truth is True:
                    self._clause_states.add((n, c, a, v))

            #..................................................................#
            # 'node(n)' - node
            #..................................................................#
            elif symbol.name == 'node':
                n = symbol.arguments[0].string
                self._map_p2l['node'][n] = guess_lit
                self._map_l2p[guess_lit] = ('node', n)
                init.add_watch(guess_lit)
                if guess_truth is True:
                    self._node_states.add(n)

            #..................................................................#
            # 'v(t, r, v)' - B-MSS state variable
            #..................................................................#
            elif symbol.name == 'v':
                t = (symbol.arguments[0].arguments[0].number,
                     symbol.arguments[0].arguments[1].number)
                r = symbol.arguments[1].string
                v = symbol.arguments[2].number
                if (r in candidates and r in reactions) or (r in inputs):
                    self.__time_model.add(t)
                    self._map_p2l['v'].setdefault(t, {})[(r, v)] = guess_lit
                    self._map_l2p[guess_lit] = ('v', t, r, v)
                    init.add_watch(guess_lit)
                    if guess_truth is True:
                        self._mss_states.add((t, r, v))

        self.unchecked_times_reg = set(self.__time_obs)
        self.unchecked_times_mss = set(self.__time_model)

    #--------------------------------------------------------------------------#
    # Propagator functions
    #--------------------------------------------------------------------------#

    def propagate(self, changes) -> None:
        """
        """
        values_mss = [self._map_l2p[lit][1:]
                      for lit in changes if self._map_l2p[lit][0] == 'v']
        self._mss_states.update(values_mss)

        values_reg = [self._map_l2p[lit][1:]
                      for lit in changes if self._map_l2p[lit][0] == 'w']
        self._reg_states.update(values_reg)

        values_node = [self._map_l2p[lit][1]
                       for lit in changes if self._map_l2p[lit][0] == 'node']
        self._node_states.update(values_node)

        values_clause = [self._map_l2p[lit][1:]
                         for lit in changes if self._map_l2p[lit][0] == 'clause']
        self._clause_states.update(values_clause)

    def propagate_verification(self, changes) -> bool:
        """
        """
        #......................................................................#
        # Retrieve the predicate associated to each assigned literal
        #......................................................................#
        readable_changes = self.translate(changes)

        #......................................................................#
        # Retrieve the time that must be checked
        #......................................................................#
        times_reg_check = set(
            l[1] for l in readable_changes
            if l[0] == 'w'
            and self.__fully_assigned_reg(l[1])
        )

        #......................................................................#
        # Check both valid regulations and spurious B-MSS
        #......................................................................#
        reg_status, reg_nogoods = self.REG_verification(
            times_reg_check,
            self._reg_states,
            self._node_states
        )

        if reg_status and len(reg_nogoods) != 0:
            return reg_status, reg_nogoods

        times_mss_check = set(
            l[1] for l in readable_changes
            if l[0] == 'v'
            and self.__fully_assigned_mss(l[1])
            and l[1][1] != 0
        )

        mss_status, mss_nogoods = self.MSS_verification(
            times_mss_check,
            self._mss_states,
            self._node_states,
            self._clause_states
        )

        return mss_status, mss_nogoods

    def undo(self, changes) -> None:
        """
        """
        values_mss = [self._map_l2p[lit][1:]
                      for lit in changes if self._map_l2p[lit][0] == 'v']
        self._mss_states.difference_update(values_mss)

        values_reg = [self._map_l2p[lit][1:]
                      for lit in changes if self._map_l2p[lit][0] == 'w']
        self._reg_states.difference_update(values_reg)

        values_node = [self._map_l2p[lit][1]
                       for lit in changes if self._map_l2p[lit][0] == 'node']
        self._node_states.difference_update(values_node)

        values_clause = [self._map_l2p[lit][1:]
                         for lit in changes if self._map_l2p[lit][0] == 'clause']
        self._clause_states.difference_update(values_clause)

        alter_reg_times = set([t for (t, _, _) in values_reg])
        self.unchecked_times_reg.update(alter_reg_times)

        alter_mss_times = set([t for (t, _, _) in values_mss])
        self.unchecked_times_mss.update(alter_mss_times)

    def check_verification(self, controler) -> bool:
        """
        """
        #......................................................................#
        # Retrieve the value assigned to each variable
        #......................................................................#
        reactions = self.__metabolic_network.get_reactions()
        inputs = self.__metabolic_network.get_inputs()

        mss = set(p[1:] for lit, p in self._map_l2p.items()
                  if p[0] == 'v'
                  and controler.assignment.is_true(lit)
                  and (p[2] in self.__node_candidates
                       or p[2] in inputs
                       )
                  )

        reg = set(p[1:] for lit, p in self._map_l2p.items()
                  if p[0] == 'w'
                  and controler.assignment.is_true(lit)
                  and p[2] in reactions
                  and p[1] in self.__time_obs
                  )

        node = set(p[1] for lit, p in self._map_l2p.items()
                   if p[0] == 'node'
                   and controler.assignment.is_true(lit)
                   )

        clause = set(p[1:] for lit, p in self._map_l2p.items()
                     if p[0] == 'clause'
                     and controler.assignment.is_true(lit)
                     )

        #......................................................................#
        # Retrieve the time that must be checked
        #......................................................................#
        times_reg_check = set(
            t for t in self.unchecked_times_reg
            if t[1] != 0
        )
        reg_status, reg_nogoods = self.REG_verification(
            times_reg_check.copy(),
            reg,
            node
        )

        self.unchecked_times_reg = times_reg_check

        if reg_status and len(reg_nogoods) != 0:
            return reg_status, reg_nogoods

        times_mss_check = set(
            t for t in self.unchecked_times_mss
            if t[1] != 0
        )

        mss_status, mss_nogoods = self.MSS_verification(
            times_mss_check.copy(),
            mss,
            node,
            clause
        )

        self.unchecked_times_mss = times_mss_check

        return mss_status, mss_nogoods

    #--------------------------------------------------------------------------#
    # Nogoods functions
    #--------------------------------------------------------------------------#

    def generate_nogoods_mss(self, times, mss, node, clause) -> list:
        """
        """
        reactions = self.__metabolic_network.get_reactions()

        #......................................................................#
        # Retrieve regulating and regulated reactions
        #......................................................................#
        to_consider_reaction = set(
            [r for (_, _, r, _) in clause if r in reactions] +
            [r for r in node if r in reactions]
        )

        #......................................................................#
        # Compute the nogoods
        #......................................................................#
        nogoods = [
            [self._map_p2l['v'][t][(r, v)] for (t, r, v) in mss
             if t == dt and r in to_consider_reaction]
            for dt in times
        ]
        nogoods = [l for l in nogoods if len(l) != 0]

        return nogoods

    def generate_nogoods_reg(self, times, reg, node):
        """
        """
        reactions = self.__metabolic_network.get_reactions()

        #......................................................................#
        # Iterate over each time. There is two cases:
        #......................................................................#
        activated_reactions = {}
        for (t, s) in times:
            activated_reactions[t] = []

            #..................................................................#
            # If s = 1: to many regulations
            #..................................................................#
            if s == 1:
                activated_reactions[t] = \
                    [self._map_p2l['w'][dt][(r, v)] for (dt, r, v) in reg
                     if t == dt
                     and v == 1
                     and r in self.__node_candidates
                     and self._map_p2l['w'][dt][(r, v)] != 1]

            #..................................................................#
            # If s = -1: not enough regulations
            #..................................................................#
            else:
                activated_reactions[t] = \
                    [self._map_p2l['w'][dt][(r, v)] for (dt, r, v) in reg
                     if t == dt
                     and v == -1
                     and r in self.__node_candidates
                     and self._map_p2l['w'][dt][(r, v)] != 1] + \
                    [self._map_p2l['node'][r] for r in node
                     if self._map_p2l['node'][r] != 1]

        #......................................................................#
        # Build the nogoods
        #......................................................................#
        nogoods = [activated_reactions[t]
                   for t, _ in times if len(activated_reactions[t]) != 0]

        return nogoods

    #--------------------------------------------------------------------------#
    # Linear verifications
    #--------------------------------------------------------------------------#

    def __fully_assigned_reg(self, time) -> bool:
        """
        """
        #......................................................................#
        # the regulations can be checked when the regulation states of each
        # reaction are defined
        #......................................................................#
        reactions = self.__metabolic_network.get_reactions()
        assigned_reg_t = set(
            [r for (t, r, _) in self._reg_states if t == time])
        return assigned_reg_t == reactions

    def __fully_assigned_mss(self, time) -> bool:
        """
        """
        #......................................................................#
        # the B-MSS can be checked when all the reactions being in the PKN have
        # been assigned to a state
        #......................................................................#
        reactions = self.__metabolic_network.get_reactions()
        candidate_node_reaction = set(
            r for r in self.__node_candidates)
        assigned_mss_t = set(r for (t, r, _) in self._mss_states if t == time)
        return assigned_mss_t == candidate_node_reaction

    def __generate_reference(self, times, obs, reg) -> tuple:
        """
        """
        reactions = self.__metabolic_network.get_reactions()
        metabolites = self.__metabolic_network.get_metabolites()

        #......................................................................#
        # Solve a rFBA for each time while forcing the flux on some reactions
        #......................................................................#
        fba_model, fba_variables = get_base_fba_lp_model(
            times,
            metabolites,
            reactions,
            self.__bounds,
            self.__metabolic_network.get_stoichiometric_coefficients(),
            self.__obj)
        fba_model = add_inhibition_constraints(
            fba_model,
            fba_variables,
            reactions,
            reg)
        fba_model = add_forced_state(
            fba_model,
            fba_variables,
            reactions,
            obs)
        _, values = solve_rfba(fba_model, fba_variables, self.__obj)

        #......................................................................#
        # Compute an upper bounds of the import reactions
        #......................................................................#
        bounds = {}
        for t in times:
            for r, m in self.__metabolic_network.get_input_reactions():
                if self.__observations.loc[t][m] == 0:
                    bounds[(t, r)] = (0.0, 0.0)
                else:
                    l_bound, u_bound = self.__bounds[(t, r)]
                    ref_u_bound = fba_variables[(t, r)].varValue
                    if ref_u_bound == 0:
                        bounds[(t, r)] = (l_bound, u_bound)
                    else:
                        bounds[(t, r)] = (l_bound, ref_u_bound)

        #......................................................................#
        # Return the objective function value of each time (to used as ref)
        #......................................................................#
        return values, bounds

    @count_and_time
    def MSS_verification(self, times, mss, node, clause) -> tuple:
        """
        """
        reactions = self.__metabolic_network.get_reactions()
        metabolites = self.__metabolic_network.get_metabolites()

        # return False, []

        #......................................................................#
        # Exit if no times must be checked
        #......................................................................#
        if len(times) == 0:
            self.nogoods = None
            return False, []

        #......................................................................#
        # Update the set of unchecked times
        #......................................................................#
        self.unchecked_times_mss.difference_update(times)

        #......................................................................#
        # Retrieves all the reactions which are used in a regulatory rule
        #......................................................................#
        used_in_clauses = set(r for (_, _, r, _) in clause if r in reactions)

        #......................................................................#
        # Retrieves current states of regulating and regulated reactions
        #......................................................................#
        mss_states = set((t, r, v) for (t, r, v) in mss
                         if t in times and (r in used_in_clauses or r in node))

        ########################################################################
        # Cache
        #......................................................................#
        cache_mss = {dt: tuple((t, r, v) for (t, r, v) in mss_states)
                     for dt in times}
        for dt, cache_dt in cache_mss.items():
            if cache_dt in self._cache_valid:
                times.remove(dt)
        if len(times) == 0:
            return False, []
        ########################################################################
        input_reactions = {
            r: m
            for r, m in self.__metabolic_network.get_input_reactions()
        }
        internal_bounds = {
            (t, r): self.__metabolic_network.get_reaction_bounds(r)
            for t in times
            for r in reactions
            if r not in input_reactions
        }
        external_bounds = {
            (t, r): self.__metabolic_network.get_reaction_bounds(r)
            if (t, m, 1) in mss else (0, 0)
            for t in times
            for r, m in input_reactions.items()
        }
        bounds = internal_bounds | external_bounds

        #......................................................................#
        # Build and solve the linear model to detect spurious B-MSS
        #......................................................................#
        fba_model, fba_variables = get_base_fba_lp_model(
            times,
            metabolites,
            reactions,
            bounds,
            self.__metabolic_network.get_stoichiometric_coefficients(),
            self.__obj,
            removeObj=True)
        fba_model, y = add_maximised_forced_state(
            fba_model,
            fba_variables,
            reactions,
            mss_states,
            bounds)
        _, values = solve_lp_model(fba_model, y)

        ########################################################################
        # Update the cache
        #......................................................................#
        for dt in values[0]:
            self._cache_valid.add(cache_mss[dt])
        ########################################################################

        #......................................................................#
        # If there is at least one spurious B-MSS, generate a nogood
        #......................................................................#
        spurious = values[1]
        if len(spurious) != 0:
            nogoods = self.generate_nogoods_mss(
                spurious,
                mss,
                node,
                clause
            )
            if len(nogoods) != 0:
                return True, nogoods

        return False, []

    @count_and_time
    def REG_verification(self, times, reg, node) -> tuple:
        """
        """
        reactions = self.__metabolic_network.get_reactions()
        metabolites = self.__metabolic_network.get_metabolites()
        #......................................................................#
        # Exit if no times must be checked
        #......................................................................#
        if len(times) == 0:
            self.nogoods = None
            return False, []

        #......................................................................#
        # Update the set of unchecked times
        #......................................................................#
        self.unchecked_times_reg.difference_update(times)

        #......................................................................#
        # Retrieve the regulatory states
        #......................................................................#
        reg_states = [(t, r, v) for (t, r, v) in reg
                      if r in node
                      and t in times
                      and r in self.__node_candidates]

        #......................................................................#
        # Generate missing reference values
        #......................................................................#
        times_no_ref = set(
            t for t in times
            if self.__objective_values.get(t, None) is None
        )
        obs_states = [
            (t, r, 1 if self.__observations.loc[t][r] > THRESHOLD else -1)
            for t in times_no_ref for r in reactions
            if r in self.__observations.loc[t]
            and self.__observations.loc[t][r] is not None
        ]
        reg_states_noref = [
            (t, r, v) for (t, r, v) in reg_states
            if t in times_no_ref
            and r in self.__node_candidates
        ]
        gen_ref_values, gen_ref_bounds = self.__generate_reference(
            times_no_ref,
            obs_states,
            reg_states_noref
        )

        bounds = self.__bounds | gen_ref_bounds

        #......................................................................#
        # Build and solve a linear model computing a rFBA
        #......................................................................#
        fba_model, fba_variables = get_base_fba_lp_model(
            times,
            metabolites,
            reactions,
            bounds,
            self.__metabolic_network.get_stoichiometric_coefficients(),
            self.__obj)
        fba_model = add_inhibition_constraints(
            fba_model,
            fba_variables,
            reactions,
            reg_states)
        _, values = solve_rfba(fba_model, fba_variables, self.__obj)

        #......................................................................#
        # For each time, compare the results of the rFBA with the reference:
        # - if reference - epsilon < rFBA < reference + epsilon: accepted
        # - else: rejected
        #......................................................................#
        unsat_times = []
        satisfiable = False
        for t in times:
            generated_ref = False
            computed_value = values[t]
            observed_value = self.__objective_values.get(t, None)
            if observed_value is None:
                observed_value = gen_ref_values[t]
                generated_ref = True

            error_range = self.__data_error if not generated_ref else 0
            is_greater = computed_value > observed_value / (1 - error_range)
            if is_greater:
                unsat_times.append((t, 1))

            is_smaller = computed_value < observed_value / (1 + error_range)
            if is_smaller:
                unsat_times.append((t, -1))

            satisfiable = satisfiable or (is_greater or is_smaller)

        #......................................................................#
        # If there is at least one rejected time, generate a nogood
        #......................................................................#
        if satisfiable:
            nogoods = self.generate_nogoods_reg(
                unsat_times,
                reg,
                node
            )
            if len(nogoods) != 0:
                return True, nogoods

        return False, []

    #--------------------------------------------------------------------------#
    # Auxiliary functions
    #--------------------------------------------------------------------------#

    def translate(self, lits) -> list:
        """
        """
        return [self._map_l2p[lit] for lit in lits]

    def nogood2str(self, nogood) -> str:
        """
        """
        predicates = []
        for lit in nogood:
            assert(lit in self._map_l2p)
            if self._map_l2p[lit][0] == 'v':
                t = self._map_l2p[lit][1]
                r = self._map_l2p[lit][2]
                v = self._map_l2p[lit][3]
                predicates.append(f'v({t}, \"{r}\", {v})')
            elif self._map_l2p[lit][0] == 'w':
                t = self._map_l2p[lit][1]
                r = self._map_l2p[lit][2]
                v = self._map_l2p[lit][3]
                predicates.append(f'w_obs({t}, \"{r}\", {v})')
            elif self._map_l2p[lit][0] == 'node':
                n = self._map_l2p[lit][1]
                predicates.append(f'node(\"{n}\")')
            elif self._map_l2p[lit][0] == 'clause':
                n = self._map_l2p[lit][1]
                c = self._map_l2p[lit][2]
                a = self._map_l2p[lit][3]
                v = self._map_l2p[lit][4]
                predicates.append(f'node(\"{n}\", {c}, \"{a}\", {v})')
        return ':-' + ', '.join(predicates) + '.'

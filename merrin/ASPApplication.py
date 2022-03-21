#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################
from lib2to3.pgen2.token import OP
import clingo
import networkx as nx
from enum import Enum

from . import ASPModel
from .Parameters import DATA_ERROR_EPSILON, MAX_OBS_BUFFER, UpdateMode
from .FBAPropagator import FBAPropagator
from .MetabolicNetworks import MetabolicNetwork

#ENUM###########################################################################


class DisplayMode(Enum):
    NONE = 'none'
    ALL_SOLUTIONS = 'all solutions'
    DEBUG = 'debug'


class ProjectionMode(Enum):
    NODE = 'node'
    COMPONENT = 'connected component'
    NETWORK = 'Boolean network'


class ResolutionMode(Enum):
    ALL = 'all'
    SUBSET_MINIMAL = 'subset_minimal'


class OptimisationMode(Enum):
    NONE = 'none'
    TIME_BUFFER = 'time buffer'


#CLASS##########################################################################


class Model(object):

    def __init__(self,
                 metabolic_network: MetabolicNetwork,
                 pkn,
                 objective,
                 n=0,
                 thread=1,
                 statistics=False,
                 resolution_mode: ResolutionMode = ResolutionMode.SUBSET_MINIMAL,
                 optimisation_mode: OptimisationMode = OptimisationMode.NONE,
                 data_error: float = DATA_ERROR_EPSILON
                 ):

        self.__metabolic_network = metabolic_network
        self.__pkn = pkn
        self.__objective = objective

        self.__results = []
        self.__status = None
        self.__asp_model = ['#program base.']
        self.__data_error = data_error

        # Initialise the controller
        options = [f'-n {n}', f'-t {thread}', '--project',
                   '--opt-strategy=usc', '--heuristic=Domain',
                   '-c bounded_nonreach=0']
        if resolution_mode == ResolutionMode.SUBSET_MINIMAL:
            options = options + ['--enum-mode=domRec', '--dom-mod=5,16']
        if optimisation_mode == OptimisationMode.TIME_BUFFER:
            options = options + ['--opt-mode=optN']
            self.__is_optimisation_problem = True
        else:
            self.__is_optimisation_problem = False
        if statistics:
            options.append("--stats")

        self.__options = options
        self.__ctl = clingo.Control(self.__options)

    def __save_assignment(self, model: clingo.Model) -> None:
        if model.optimality_proven or not self.__is_optimisation_problem:
            show_terms = model.symbols(shown=True)

            # Term of the form: clause(N, C, A, V)
            terms = {n: {} for n in self.__pkn.nodes}
            for term in show_terms:
                N = term.arguments[0].string
                C = term.arguments[1].number
                A = term.arguments[2].string
                V = term.arguments[3].number
                terms[N].setdefault(C, []).append((A, V))

            result = {}
            for n in terms.keys():
                rule = []
                for c in terms[n].keys():
                    terms[n][c].sort()
                    clause = [('!' if v == -1 else '') +
                              m for (m, v) in terms[n][c]]
                    clause = ' & '.join(clause)
                    rule.append(clause)
                if len(rule) > 1:
                    rule = [f'({c})' for c in rule]
                rule = ' || '.join(rule)
                result[n] = rule

            self.__results.append(result)

    def build(self, observations, update_mode=UpdateMode.SYNCHRONOUS, observation_buffer: int = MAX_OBS_BUFFER):
        """
        """
        # Add the linear propagator
        self.__propagator = FBAPropagator(
            self.__metabolic_network,
            self.__pkn,
            observations,
            self.__objective,
            data_error = self.__data_error)
        self.__ctl.register_propagator(self.__propagator)

        # Data
        self.custom(ASPModel.asp_metabolic_network(self.__metabolic_network))
        self.custom(ASPModel.asp_pkn(self.__pkn))

        regulators = set(n for n in self.__pkn.nodes
                         if n not in self.__metabolic_network.get_reactions()
                         and n not in self.__metabolic_network.get_inputs()
                         and n not in self.__metabolic_network.get_outputs())
        self.custom(ASPModel.asp_fixed_nodes(regulators))

        self.custom(ASPModel.asp_observations(
            observations, max_obs_buffer=observation_buffer))

        # Inferring Model
        self.custom(ASPModel.inferring_model)

        # Options
        self.custom(ASPModel.update_mode(update_mode))

        # Optimisation:
        if '--opt-mode=optN' in self.__options:
            self.custom(ASPModel.optimisation_time_buffer())

    def custom(self, asp: str) -> None:
        self.__ctl.add("base", [], asp)
        self.__asp_model.append(asp)

    def forced_regulation(self, node) -> None:
        if type(node) is list:
            for n in node:
                self.custom(f'node("{n}").')
        else:
            self.custom(f'node("{node}").')

    def ground(self) -> None:
        self.__ctl.ground([('base', [])])

    def solve(self, timelimit: int = None) -> clingo.SolveResult:
        self.custom(ASPModel.show_projection_network())

        self.ground()

        self.__results = []

        # Solve the model        
        with self.__ctl.solve(on_model=self.__save_assignment, async_=True) as handle:
            handle.wait(timelimit)
            handle.cancel()
            self.__status = handle.get()

        return self.__status

    def statistics(self):
        stats  = {}
        for meth in ["MSS_verification", "REG_verification"]:
            stats[meth] = {"nb_calls": 0, "total_duration": 0}
            for prop in self.__propagator._fluxprops:
                m = getattr(prop, f"__stats_{meth}")
                stats[meth]["nb_calls"] += m["calls"]
                stats[meth]["total_duration"] += m["total_duration"]
        stats["clingo"] = self.__ctl.statistics
        return stats

    def __solve_projection_node(self) -> dict:
        """
        """
        asp_model = '\n'.join(self.__asp_model)
        results = {}
        for node in self.__pkn.nodes:
            if len(list(self.__pkn.predecessors(node))) == 0:
                results[node] = []
                continue

            self.__ctl = clingo.Control(self.__options)
            self.__propagator.reset()
            self.__ctl.register_propagator(self.__propagator)

            self.__ctl.add("base", [], asp_model)
            self.__ctl.add("base", [], ASPModel.show_projection_node(node))
            self.ground()

            self.__results = []

            # Solve the model
            self.__status = self.__ctl.solve(on_model=self.__save_assignment)
            if len(self.__propagator.generated_nogoods) != 0:
                asp_model += '\n' + \
                    '\n'.join(self.__propagator.generated_nogoods)
                self.__propagator.reset_memory()

            # TODO: traiter le cas UNSAT
            # Normalement : si un est UNSAT, tous les autres le sont.

            for answer_set in self.get_results():
                results.setdefault(node, []).append(answer_set[node])
        return results

    def __solve_projection_component(self) -> dict:
        """
        """
        asp_model = '\n'.join(self.__asp_model)
        results = {}
        id = 1
        for component in nx.weakly_connected_components(self.__pkn):
            self.__ctl = clingo.Control(self.__options)

            self.__propagator.reset()
            self.__ctl.register_propagator(self.__propagator)

            self.__ctl.add("base", [], asp_model)
            self.__ctl.add(
                "base", [], ASPModel.show_projection_component(component))
            self.ground()

            self.__results = []

            # Solve the model
            self.__status = self.__ctl.solve(on_model=self.__save_assignment)
            if len(self.__propagator.generated_nogoods) != 0:
                asp_model += '\n' + \
                    '\n'.join(self.__propagator.generated_nogoods)
                self.__propagator.reset_memory()

            # TODO: traiter le cas UNSAT
            # Normalement : si un est UNSAT, tous les autres le sont.

            for answer_set in self.get_results():
                results.setdefault(id, []).append(
                    {k: v for k, v in answer_set.items() if k in component})
            id += 1
        return results

    def solve_projection(self, projection_mode: ProjectionMode) -> dict:
        """
        """
        if projection_mode == ProjectionMode.NETWORK:
            self.solve()
            results = self.get_results()
        elif projection_mode == ProjectionMode.COMPONENT:
            results = self.__solve_projection_component()
        elif projection_mode == ProjectionMode.NODE:
            results = self.__solve_projection_node()
        return results

    def get_results(self) -> list:
        return self.__results.copy()

    def export_model(self, nogoods=False) -> str:
        asp_model = '\n'.join(self.__asp_model)
        if nogoods:
            if len(self.__propagator.generated_nogoods) != 0:
                nogoods = '\n' + '\n'.join(self.__propagator.generated_nogoods)
                asp_model += nogoods
        return asp_model

    def is_satisfiable(self) -> bool:
        return self.__status.satisfiable

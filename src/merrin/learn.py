# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from typing import Literal, Iterable
from enum import Enum

from pandas import DataFrame
from clingo import Control, Model

from merrinasp.theory.language import THEORY_LANGUAGE, rewrite
from merrinasp.theory.propagator import LpPropagator

from merrin.datastructure import MetabolicNetwork, Observation
from merrin.asp import (
    instantiate_networks,
    instantiate_observations,
    instantiate_parameters,
    ASP_MODEL_LEARN
)

# ==============================================================================
# Class: MerrinLearner
# ==============================================================================
# Object used to infer metabolic regulatory rules from timeseries data
class MerrinLearner:

    def __init__(self: MerrinLearner, sbml: str, objective: str,
                 pkn: Iterable[tuple[str, int, str]], max_clause: int = 20,
                 mn_bounds: Iterable[tuple[str, float, float]] | None = None) \
            -> None:
        # ----------------------------------------------------------------------
        # Load metabolic network
        # ----------------------------------------------------------------------
        self.__mn: MetabolicNetwork = MetabolicNetwork.read_sbml(sbml)
        if mn_bounds is not None:
            for r, lb, ub in mn_bounds:
                self.__mn.set_bound(r, lb, ub)
        self.__mn = self.__mn.to_irreversible()
        self.__renamed_reactions: dict[str, tuple[str, str]] = {
            r_: (r, f'[{r} < 0]') if r_ == rr else (r, f'[{r} > 0]')
            for r, (rr, rf) in \
                self.__mn.previously_reversible_reactions().items()
            for r_ in (rr, rf)
        }
        node_to_rename: dict[str, tuple[str, str]] = \
            self.__mn.previously_reversible_reactions()
        # ~ Transfert the MN irreversibility conversion to the objective
        self.__objective: str = objective if objective not in node_to_rename \
                                    else node_to_rename[objective][1]
        # ~ Transfert the MN irreversibility conversion to the PKN
        self.__pkn: set[tuple[str, int, str]] = set()
        for u, s, v in pkn:
            u_l: Iterable[str] = \
                [u] if u not in node_to_rename else node_to_rename[u]
            v_l: Iterable[str] = \
                [v] if v not in node_to_rename else node_to_rename[v]
            for u_ in u_l:
                for v_ in v_l:
                    self.__pkn.add((u_, s, v_))
        # ----------------------------------------------------------------------
        # Instantiate ASP constraints
        # ----------------------------------------------------------------------
        self.__constraints: list[str] = []
        self.__constraints.extend(instantiate_networks(
            self.__mn,
            objective,
            pkn,
            max_clause=max_clause
        ))
        self.__projection_constraints: list[str] = []
        self.__observations_constraints: list[str] = []
        # ----------------------------------------------------------------------
        # Set default solving parameters
        # ----------------------------------------------------------------------
        self.__propagator: LpPropagator | None = None
        self.__projection: MerrinLearner.Projection = \
            MerrinLearner.Projection.NETWORK
        self.__optimisation: MerrinLearner.Optimisation = \
            MerrinLearner.Optimisation.SUBSETMIN

    # ==========================================================================
    # Parameters
    # ==========================================================================
    class Projection(Enum):
        NETWORK = 0
        NODE = 1

    class Optimisation(Enum):
        ALL = 1
        SUBSETMIN = 2

    # ==========================================================================
    # Solving
    # ==========================================================================
    def learn(self: MerrinLearner, observations: Iterable[Observation],
              nbsol: int = 0, max_gap: int = 10, max_error: float = 0.1,
              lp_epsilon: float = 10**-5, timelimit: int = -1,
              lpsolver: Literal['glpk', 'gurobi'] = 'glpk') \
            -> DataFrame:
        # ----------------------------------------------------------------------
        # Instantiate missing ASP constraints
        # ----------------------------------------------------------------------
        self.__observations_constraints.clear()
        self.__projection_constraints.clear()

        # ~ Observation constraints
        self.__observations_constraints.extend(
            instantiate_observations(
                self.__mn, self.__objective, observations, epsilon=10**-9
            )
        )
        # ~ Parameters constraints
        self.__observations_constraints.extend(
            instantiate_parameters(max_gap, max_error, lp_epsilon)
        )
        # ----------------------------------------------------------------------
        # Launch the solving process
        # ----------------------------------------------------------------------
        if self.__projection == MerrinLearner.Projection.NETWORK:
            return self.__learn_rn(nbsol, timelimit, lpsolver)
        if self.__projection == MerrinLearner.Projection.NODE:
            return self.__learn_nodes(nbsol, timelimit, lpsolver)
        assert False

    def __learn_rn(self: MerrinLearner, nbsol: int = 0, timelimit: int = -1,
                   lpsolver: Literal['glpk', 'gurobi'] = 'glpk') \
            -> DataFrame:
        # ----------------------------------------------------------------------
        # Initialise the ASP solver `clingo`
        # ----------------------------------------------------------------------
        # ~ Solving options
        options: list[str] = self.__get_options(nbsol)
        # Initialize `clingo` controller
        ctl = Control(options)
        self.__propagator = LpPropagator(lpsolver)
        self.__propagator.lazy(False)
        self.__propagator.show_lpassignment(False)
        self.__propagator.strict_forall_check(False)
        ctl.register_propagator(self.__propagator)  # type: ignore
        # ----------------------------------------------------------------------
        # Build and solve the ASP program
        # ----------------------------------------------------------------------
        # ~ Load the theory language
        ctl.add("base", [], THEORY_LANGUAGE)
        # ~ Add all constraints
        ctl.add("base", [], '\n'.join(self.__constraints +
                                      self.__observations_constraints +
                                      self.__projection_constraints))
        # ~ Load the model
        rewrite(ctl, [ASP_MODEL_LEARN])
        # ~ Add Show constraint
        ctl.add("base", [], '\n'.join([
            '#show.',
            '#show clause/4.'
        ]))
        # ~ Ground the ASP program
        ctl.ground([('base', [])])
        # ~ Solve
        nodes: set[str] = {n for _, _, n in self.__pkn}
        results: DataFrame = DataFrame(columns=list(sorted(nodes)))
        if timelimit == -1:
            ctl.solve(on_model=lambda m: self.__on_model_rn(results, m))
        else:
            with ctl.solve(on_model=lambda m: self.__on_model_rn(results, m),
                           async_=True) as handle:  # type: ignore
                handle.wait(timelimit)
                handle.cancel()
                handle.get()
        return results

    def __learn_nodes(self: MerrinLearner, nbsol: int = 0, timelimit: int = -1,
                      lpsolver: Literal['glpk', 'gurobi'] = 'glpk') \
            -> DataFrame:
        nodes: set[str] = {n for _, _, n in self.__pkn}
        results: DataFrame = DataFrame(columns=list(sorted(nodes)))
        options: list[str] = self.__get_options(nbsol)
        for n in sorted(nodes):
            # Initialize `clingo` controller
            ctl = Control(options)
            self.__propagator = LpPropagator(lpsolver)
            self.__propagator.lazy(False)
            self.__propagator.show_lpassignment(False)
            self.__propagator.strict_forall_check(False)
            ctl.register_propagator(self.__propagator)  # type: ignore
            # Solve for the node `n`
            n_rules = self.__learn_node(ctl, n, timelimit=timelimit)
            results[n] = [';'.join(sorted(n_rules))]
        return results

    def __learn_node(self: MerrinLearner, ctl: Control, n: str,
                     timelimit: int = -1) -> list[str]:
        # ~ Load the theory language
        ctl.add("base", [], THEORY_LANGUAGE)
        # ~ Add all constraints
        ctl.add("base", [], '\n'.join(self.__constraints +
                                      self.__observations_constraints +
                                      self.__projection_constraints))
        # ~ Load the model
        rewrite(ctl, [ASP_MODEL_LEARN])
        # ~ Add Show constraint
        ctl.add("base", [], '\n'.join([
            '#show.',
            f'#show clause(N,C,A,V): clause(N,C,A,V), N="{n}".'
        ]))
        # ~ Ground the ASP program
        ctl.ground([('base', [])])
        # ~ Solve
        results: list[str] = []
        if timelimit == -1:
            ctl.solve(on_model=lambda m: self.__on_model_node(results, m))
        else:
            with ctl.solve(on_model=lambda m: self.__on_model_node(results, m),
                           async_=True) as handle:  # type: ignore
                handle.wait(timelimit)
                handle.cancel()
                handle.get()
        return results

    # --------------------------------------------------------------------------
    # Auxiliary ASP functions
    # --------------------------------------------------------------------------
    def __on_model_rn(self: MerrinLearner, results: DataFrame,
                      model: Model) -> None:
        if not model.optimality_proven:
            return
        nodes: set[str] = {n for _, _, n in self.__pkn}
        rules: dict[str, str] = self.__parse_clauses(model)
        results.loc[len(results.index)] = \
            [rules.get(n, '1') for n in sorted(nodes)]

    def __on_model_node(self: MerrinLearner, results: list[str],
                        model: Model) -> None:
        if not model.optimality_proven:
            return
        rules: dict[str, str] = self.__parse_clauses(model)
        if len(rules) == 0:
            results.append('1')
        else:
            for _, v in rules.items():
                results.append(v)

    # ==========================================================================
    # Clingo output parsing
    # ==========================================================================
    def __parse_clauses(self: MerrinLearner, model: Model) -> dict[str, str]:
                # ~ Extract clauses from ASP solutions
        clauses: dict[str, dict[int, list[str]]] = {}
        for atom in model.symbols(shown=True):
            assert atom.name == 'clause' and len(atom.arguments) == 4
            n: str = atom.arguments[0].string
            c: int = atom.arguments[1].number
            a: str = atom.arguments[2].string
            if a in self.__renamed_reactions:
                a = self.__renamed_reactions[a][1]
            v: int = atom.arguments[3].number
            clauses.setdefault(n, {}).setdefault(c, []).append(
                a if v == 1 else f'!{a}'
            )
        # ~ Generate rule from the clauses
        rules: dict[str, str] = {}
        for n, n_clauses in clauses.items():
            rule_list: list[str] = []
            for c in sorted(n_clauses):
                assert len(n_clauses[c]) > 0
                clause: str = ' & '.join(n_clauses[c])
                if len(n_clauses[c]) > 1:
                    clause = f'({clause})'
                rule_list.append(clause)
            rule: str = ' | '.join(rule_list)
            if len(rule_list) > 1:
                rule = f'({rule})'
            if n in self.__renamed_reactions:
                n = self.__renamed_reactions[n][0]
            rules[n] = rule
        return rules

    # ==========================================================================
    # Setters / Getters
    # ==========================================================================
    # --------------------------------------------------------------------------
    # Setters
    # --------------------------------------------------------------------------
    def set_projection(self: MerrinLearner, mode: MerrinLearner.Projection) \
            -> None:
        self.__projection = mode

    def set_optimisation(self: MerrinLearner,
                         mode: MerrinLearner.Optimisation) -> None:
        self.__optimisation = mode

    # --------------------------------------------------------------------------
    # Getters
    # --------------------------------------------------------------------------
    def __get_options(self: MerrinLearner, nbsol: int = 0) -> list[str]:
        options: list[str] = [
            f'-n {nbsol}', '-t 1', '--project',
            '--opt-mode=optN', '--opt-strategy=usc',
            '--heuristic=Domain', '-c bounded_nonreach=0'
        ]
        if self.__optimisation == MerrinLearner.Optimisation.SUBSETMIN:
            options += ['--enum-mode=domRec', '--dom-mod=5,16']
        return options

    def statistics(self: MerrinLearner) -> None:
        raise NotImplementedError()

    def get_model(self: MerrinLearner) -> tuple[str, str]:
        cmd: str = ''  # TODO: auto compute the bash command
        model: str = '\n'.join(
            self.__constraints +
            self.__observations_constraints
        )
        with open(ASP_MODEL_LEARN, 'r', encoding='utf-8') as file:
            for line in file.readlines():
                model += '\n' + line.strip()
        return (cmd, model)

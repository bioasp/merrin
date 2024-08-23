# ==============================================================================
# Imports
# ==============================================================================
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Literal, Iterable, Callable

from clingo import Control, Model

from merrinasp.theory.language import THEORY_LANGUAGE, rewrite
from merrinasp.theory.propagator import LpPropagator

from merrin.datastructure import MetabolicNetwork, Observation
from merrin.asp import (
    instantiate_mn,
    instantiate_pkn,
    instantiate_observations,
    instantiate_parameters,
    instantiate_trace,
    ASP_MODEL_LEARN,
    ASP_MODEL_LEARN_FROM_TRACE
)


# ==============================================================================
# Class: MerrinLearner
# ==============================================================================
# Object used to infer metabolic regulatory rules from timeseries data
class MerrinLearner:

    @dataclass
    class __Config:
        max_clause: int = 20
        max_gap: int = 10
        max_error: float = 0.1
        lp_epsilon: float = 10**-5
        timelimit: int = -1
        lp_solver: Literal['glpk', 'gurobi'] = 'glpk'

    # ==========================================================================
    # Initialisation
    # ==========================================================================
    def __init__(self: MerrinLearner) -> None:
        # ~ Internal state
        self.__instantiated: bool = False
        self.__propagator: LpPropagator | None = None
        # ~ ASP
        self.__constraints: list[str] = []
        self.__constraints_extended: list[str] = []
        # ~ Cache
        self.__pkn: set[tuple[str, int, str]] = set()
        self.__renamed_reactions: dict[str, tuple[str, str]] = {}

    # ==========================================================================
    # Loaders
    # ==========================================================================
    def load_instance(self: MerrinLearner, mn: MetabolicNetwork,
                      objective: str, pkn: Iterable[tuple[str, int, str]],
                      observations: Iterable[Observation],
                      epsilon: float = 10**-9) -> None:
        # ~ Pre-process Metabolic Network
        mn.to_irreversible()
        self.__renamed_reactions = {
            r_: (r, f'[{r} < 0]') if r_ == rr else (r, f'[{r} > 0]')
            for r, (rf, rr) in mn.previously_reversible_reactions().items()
            for r_ in (rf, rr)
        }
        assert mn.irreversible
        assert objective in mn.reactions()
        # ~ Save PKN
        self.__pkn.clear()
        for u, s, v in pkn:
            for u_ in mn.previously_reversible_reactions().get(u, (u,)):
                for v_ in mn.previously_reversible_reactions().get(v, (v,)):
                    self.__pkn.add((u_, s, v_))

        # ~ Build ASP constraints
        self.__build_asp(mn, objective, observations, epsilon)
        # ~ Update the instantiated status
        self.__instantiated = True

    def __build_asp(self: MerrinLearner, mn: MetabolicNetwork, objective: str,
                    observations, epsilon) -> None:
        self.__constraints.clear()
        self.__constraints.extend(instantiate_mn(mn, objective))
        self.__constraints.extend(
            instantiate_observations(mn, objective, observations, epsilon)
        )

    def __build_asp_pkn(self: MerrinLearner,
                        config: MerrinLearner.__Config) -> None:
        self.__constraints_extended.clear()
        self.__constraints_extended.extend(
            instantiate_parameters(config.max_gap,
                                   config.max_error,
                                   config.lp_epsilon)
        )
        self.__constraints_extended.extend(
            instantiate_pkn(self.__pkn, config.max_clause)
        )

    def __read_config(self: MerrinLearner,
                      kwargs: dict[str, Any]) -> MerrinLearner.__Config:
        # ~ Parse arguments
        max_clause: int = 20
        max_gap: int = kwargs.get('max_gap', 10)
        max_error: float = kwargs.get('max_error', 0.1)
        lp_epsilon: float = kwargs.get('lp_epsilon', 10**-5)
        timelimit: int = kwargs.get('timelimit', -1)
        lp_solver: Literal['glpk', 'gurobi'] = 'glpk'
        if kwargs.get('lp_solver', 'glpk') in ['glpk', 'gurobi']:
            lp_solver = kwargs.get('lp_solver', 'glpk')
        # ~ Build dataclass
        return MerrinLearner.__Config(max_clause, max_gap, max_error,
                                      lp_epsilon, timelimit, lp_solver)

    def __init_clingo(self: MerrinLearner, nbsol: int, subsetmin: bool,
                      lp_solver: Literal['glpk', 'gurobi']) -> Control:
        # ~ Solving options
        options: list[str] = self.__get_options(nbsol, subsetmin)
        # ~ Propagator
        self.__propagator = LpPropagator(lp_solver)
        self.__propagator.lazy(False)
        self.__propagator.show_lpassignment(False)
        self.__propagator.strict_forall_check(False)
        # ~ `clingo` Controller
        ctl = Control(options)
        ctl.register_propagator(self.__propagator)  # type: ignore
        # ~ Load the theory language and ASP model
        ctl.add("base", [], THEORY_LANGUAGE)
        rewrite(ctl, [ASP_MODEL_LEARN])
        ctl.add("base", [], '\n'.join(self.__constraints +
                                      self.__constraints_extended))
        return ctl

    def __solve_asp(self: MerrinLearner, ctl: Control, timelimit: int,
                    on_model: Callable[[Model], None]) -> None:
        if timelimit == -1:
            ctl.solve(on_model=on_model)
            return

        with ctl.solve(on_model=on_model, async_=True) as handle:
            handle.wait(timelimit)
            handle.cancel()
            handle.get()

    def __prohibit_model(self: MerrinLearner,
                         signatures: list[tuple[str, int]],
                         model: Model) -> list[tuple[str, bool]]:
        to_prohibit: list[tuple[str, bool]] = []
        for symbol in model.symbols(atoms=True):
            if (symbol.name, len(symbol.arguments)) in signatures:
                to_prohibit.append((str(symbol), False))
        assert self.__propagator is not None
        self.__propagator.add_clause(to_prohibit)
        return to_prohibit

    # ==========================================================================
    # Rules projection mode
    # ==========================================================================

    # --------------------------------------------------------------------------
    # Learn: `Network` Projection
    # --------------------------------------------------------------------------
    def learn(self: MerrinLearner, nbsol: int = 0, subsetmin: bool = False,
              display: bool = False, **kwargs) -> list[list[tuple[str, str]]]:
        assert self.__instantiated
        # ~ Get Config
        config: MerrinLearner.__Config = self.__read_config(kwargs)
        # ~ Build Parameters and PKN ASP
        self.__build_asp_pkn(config)
        # ~ Initialise the ASP solver `clingo`
        ctl: Control = self.__init_clingo(nbsol, subsetmin, config.lp_solver)
        ctl.add("base", [], '\n'.join([
            '#show.',
            '#show clause/4.'
        ]))
        # ~ Ground the ASP program
        ctl.ground([('base', [])])
        # ~ Print CSV header
        if display:
            nodes: set[str] = {n for _, _, n in self.__pkn}
            columns: list[str] = [
                self.__renamed_reactions.get(n, (n, n))[1]
                for n in sorted(nodes)
            ]
            print(','.join(columns), flush=True)
        # ~ Solve
        results: list[list[tuple[str, str]]] = []
        self.__solve_asp(
            ctl, config.timelimit,
            lambda m: self.__on_model_network(results, m, display)
        )
        return results

    def __on_model_network(self: MerrinLearner,
                           results: list[list[tuple[str, str]]],
                           model: Model, display: bool = False) -> None:
        if not model.optimality_proven:
            return
        nodes: set[str] = {n for _, _, n in self.__pkn}
        rules: dict[str, str] = self.__parse_clauses(model)
        results.append(sorted(rules.items()))
        if display:
            print(','.join([rules.get(n, '1') for n in sorted(nodes)]),
                  flush=True)

    # --------------------------------------------------------------------------
    # Learn: `Node` Projection
    # --------------------------------------------------------------------------
    def learn_per_node(self: MerrinLearner, nbsol: int = 0,
                       subsetmin: bool = False, display: bool = False,
                       **kwargs) -> dict[str, list[str]]:
        assert self.__instantiated
        # ~ Get Config
        config: MerrinLearner.__Config = self.__read_config(kwargs)
        # ~ Build Parameters and PKN ASP
        self.__build_asp_pkn(config)
        # ~ Initialise the ASP solver `clingo`
        nodes: set[str] = {n for _, _, n in self.__pkn}
        # ~ Print CSV header
        if display:
            columns: list[str] = [
                self.__renamed_reactions.get(n, (n, n))[1]
                for n in sorted(nodes)
            ]
            print(','.join(columns), flush=True)
        ctl: Control = self.__init_clingo(nbsol, subsetmin,
                                            config.lp_solver)
        # ~ Add show
        ctl.add("base", [], '\n'.join([
                '1 { show(N): in(_,N,_) } 1.',
                '#show.',
                '#show show/1.',
                '#show clause(N,C,A,V): clause(N,C,A,V), show(N).',
            ]))
        # ~ Ground the ASP program
        ctl.ground([('base', [])])
        # ~ Solve
        results: dict[str, list[str]] = {}
        self.__solve_asp(ctl, config.timelimit,
                         lambda m: self.__on_model_node(results, m))

        if display:
            results_str: list[str] = [
                ';'.join(sorted(results[n])) for n in sorted(results.keys())
            ]
            print(','.join(results_str), flush=True)

        return results

    def __on_model_node(self: MerrinLearner,
                        results: dict[str, list[str]],
                        model: Model) -> None:
        if not model.optimality_proven:
            return
        if len(model.symbols(shown=True)) == 1:
            for atom in model.symbols(shown=True):
                assert atom.name == 'show' and len(atom.arguments) == 1
                results.setdefault(atom.arguments[0].string, []).append('1')
            return

        rules: dict[str, str] = self.__parse_clauses(model)
        assert len(rules) == 1
        for node, rule in rules.items():
            results.setdefault(node, []).append(rule)

    # --------------------------------------------------------------------------
    # Learn: `Trace` Projection
    # --------------------------------------------------------------------------
    def learn_per_trace(self: MerrinLearner, nbsol: int = 0,
                        subsetmin: bool = False, display: bool = False,
                        **kwargs) -> list[dict[str, list[str]]]:
        assert self.__instantiated
        # ~ Get Config
        config: MerrinLearner.__Config = self.__read_config(kwargs)
        # ~ Build Parameters and PKN ASP
        self.__build_asp_pkn(config)
        # ~ Initialise the ASP solver `clingo`
        ctl: Control = self.__init_clingo(nbsol, subsetmin, config.lp_solver)
        ctl.add("base", [], '\n'.join([
            'inp(T2,N,V) :- succ(T1,T2), read(T1,N,V), in(N,_,_).',
            'out(T,N,V) :- succ(_,T), x(T,N,V), in(_,N,_).',
            '#show.',
            '#show clause/4.'
        ]))
        # ~ Ground the ASP program
        ctl.ground([('base', [])])
        # ~ Print CSV header
        if display:
            nodes: set[str] = {n for _, _, n in self.__pkn}
            columns: list[str] = [
                self.__renamed_reactions.get(n, (n, n))[1]
                for n in sorted(nodes)
            ]
            print(','.join(columns), flush=True)
        # ~ Solve
        results: list[dict[str, list[str]]] = []
        self.__solve_asp(ctl, config.timelimit,
                         lambda m: self.__on_model_trace(results, m,
                                                         display, subsetmin))
        return results

    def __on_model_trace(self: MerrinLearner,
                         results: list[dict[str, list[str]]],
                         model: Model, display: bool = False,
                         subsetmin: bool = False) -> None:
        if not model.optimality_proven:
            return
        # ~ Prohibit the trace
        self.__prohibit_model([('inp', 3), ('out', 3)], model)
        if subsetmin:
            self.__prohibit_strict_supersets(model)
        # ~ Parse the trace
        trace: dict[tuple[str, int],
                    tuple[list[tuple[str, bool]],
                          list[tuple[str, bool]]]] = self.__parse_trace(model)
        # ~ Learn the BN compatible with the trace
        rules: dict[str, list[str]] = self.__learn_from_trace(trace, subsetmin)
        results.append(rules)
        # ~ Display the result
        if display:
            rules_str: list[str] = [
                ';'.join(sorted(rules[n])) for n in sorted(rules.keys())
            ]
            print(','.join(rules_str), flush=True)

    def __prohibit_strict_supersets(self: MerrinLearner, model: Model) -> None:
        nodes_true: list[tuple[str, bool]] = []
        for symbol in model.symbols(atoms=True):
            if symbol.name == 'node' and len(symbol.arguments) == 1:
                nodes_true.append((str(symbol), False))
        assert self.__propagator is not None
        for symbol in model.symbols(atoms=True, complement=True):
            if symbol.name == 'node' and len(symbol.arguments) == 1:
                self.__propagator.add_clause(
                    [(str(symbol), False)] + nodes_true
                )

    def __learn_from_trace(self: MerrinLearner,
                           trace: dict[tuple[str, int],
                                       tuple[list[tuple[str, bool]],
                                             list[tuple[str, bool]]]],
                           subsetmin: bool = False) -> dict[str, list[str]]:
        # ~ Build the ASP program
        options: list[str] = self.__get_options(0, subsetmin)
        ctl = Control(options)
        ctl.load(ASP_MODEL_LEARN_FROM_TRACE)
        ctl.add("base", [], '\n'.join(self.__constraints_extended))
        ctl.add("base", [], '\n'.join(instantiate_trace(trace)))
        # ~ For each node in the PKN with input interactions
        ctl.add("base", [], '\n'.join([
                '1 { show(N): in(_,N,_) } 1.',
                '#show.',
                '#show show/1.',
                '#show clause(N,C,A,V): clause(N,C,A,V), show(N).']))
        #   | Ground the ASP program
        ctl.ground([('base', [])])
        #   | Solve the ASP program
        results: dict[str, list[str]] = {}
        self.__solve_asp(ctl, -1,
                         lambda m: self.__on_model_trace_learn(results, m))
        # ~ Return the results
        return results

    def __on_model_trace_learn(self: MerrinLearner,
                               results: dict[str, list[str]],
                               model: Model) -> None:
        if len(model.symbols(shown=True)) == 1:
            for atom in model.symbols(shown=True):
                assert atom.name == 'show' and len(atom.arguments) == 1
                results.setdefault(atom.arguments[0].string, []).append('1')
            return

        rules: dict[str, str] = self.__parse_clauses(model)
        assert len(rules) == 1
        for node, rule in rules.items():
            results.setdefault(node, []).append(rule)

    # ==========================================================================
    # Clingo output parsing
    # ==========================================================================
    def __parse_clauses(self: MerrinLearner, model: Model) -> dict[str, str]:
        # ~ Extract clauses from ASP solutions
        clauses: dict[str, dict[int, list[str]]] = {}
        for atom in model.symbols(shown=True):
            if atom.name == 'clause':
                assert len(atom.arguments) == 4
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
            rules[n] = rule
        return rules

    def __parse_trace(self: MerrinLearner, model: Model) \
            -> dict[tuple[str, int],
                    tuple[list[tuple[str, bool]], list[tuple[str, bool]]]]:
        # ~ Extract clauses from ASP solutions
        trace: dict[tuple[str, int],
                    tuple[list[tuple[str, bool]], list[tuple[str, bool]]]] = {}
        for symbol in model.symbols(atoms=True):
            if symbol.name in ['inp', 'out']:
                index: int = 0 if symbol.name == 'inp' else 1
                time: tuple[str, int] = (
                    symbol.arguments[0].arguments[0].string,
                    symbol.arguments[0].arguments[1].number
                )
                node: str =symbol.arguments[1].string
                value: bool = symbol.arguments[2].number == 1
                trace.setdefault(time, ([], []))[index].append((node, value))
        return trace

    # ==========================================================================
    # Setters / Getters
    # ==========================================================================
    # --------------------------------------------------------------------------
    # Getters
    # --------------------------------------------------------------------------
    def __get_options(self: MerrinLearner, nbsol: int = 0,
                      subsetmin: bool = False) -> list[str]:
        options: list[str] = [
            f'-n {nbsol}', '-t 1', '--project',
            '--opt-mode=optN', '--opt-strategy=usc'
        ]
        if subsetmin:
            options += ['--enum-mode=domRec', '--dom-mod=5,16',
                        '--heuristic=Domain', '-c bounded_nonreach=0']
        return options

    def get_model(self: MerrinLearner) -> str:
        model: str = '\n'.join(
            self.__constraints +
            self.__constraints_extended
        )
        with open(ASP_MODEL_LEARN, 'r', encoding='utf-8') as file:
            for line in file.readlines():
                model += '\n' + line.strip()
        return model

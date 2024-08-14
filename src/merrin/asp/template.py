
# ==============================================================================
# Import
# ==============================================================================
from typing import Literal
from math import factorial


# ==============================================================================
# ASP Merrin instance template
# ==============================================================================

# ------------------------------------------------------------------------------
# Pretty printing
# ------------------------------------------------------------------------------
def header(title: str, sep: str) -> str:
    return '% ' + (sep * 78) + '\n' + f'% {title}\n' + '% ' + (sep * 78)


def comment(content: str) -> str:
    return f'% {content}'


# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
class Parameter:
    # --------------------------------------------------------------------------
    # Template ASP
    # --------------------------------------------------------------------------
    @classmethod
    def lp_epsilon(cls, epsilon: float) -> str:
        return f'epsilon("{epsilon}").'

    @classmethod
    def max_gap(cls, max_gap: int) -> str:
        return f'maxGap({max_gap}).'

    @classmethod
    def max_error(cls, max_error: float) -> str:
        lb: float = 1 / (1 + max_error)
        ub: float = 1 / (1 - max_error)
        return f'lb("{lb}"). ub("{ub}").'


# ------------------------------------------------------------------------------
# Metabolic network description
# ------------------------------------------------------------------------------
class MetabolicNetwork:
    # --------------------------------------------------------------------------
    # Template ASP
    # --------------------------------------------------------------------------
    @classmethod
    def metabolite_external(cls, m: str) -> str:
        return f'ext("{m}").'

    @classmethod
    def reactant(cls, r: str, m: str, s: float) -> str:
        return f'reactant("{m}", "{r}", "{s}").'

    @classmethod
    def product(cls, r: str, m: str, s: float) -> str:
        return f'product("{m}", "{r}", "{s}").'

    @classmethod
    def bounds(cls, r: str, lb: float, ub: float) -> str:
        return f'bound("{r}", "{lb}", "{ub}").'

    @classmethod
    def fixed_flux(cls, r: str) -> str:
        return f'fixedFlux("{r}").'

    @classmethod
    def reversible(cls, rf: str, rr: str) -> str:
        return f'rev("{rf}", "{rr}").'

    @classmethod
    def objective(cls, r: str) -> str:
        return f'objective("{r}").'
    
    @classmethod
    def gene(cls, g: str) -> str:
        return f'gene("{g}").'


# ------------------------------------------------------------------------------
# Prior Knowledge Network (PKN) description
# ------------------------------------------------------------------------------
class PKN:
    # --------------------------------------------------------------------------
    # Template ASP
    # --------------------------------------------------------------------------
    @classmethod
    def node(cls, n: str, fixed: bool = False) -> str:
        if fixed:
            return f'node("{n}").'
        return f'{{ node("{n}") }}.'

    @classmethod
    def edge(cls, u: str, v: str, s: int) -> str:
        if s == 0:
            return f'in("{u}", "{v}", (-1;1)).'
        return f'in("{u}", "{v}", {s}).'

    @classmethod
    def max_clause(cls, n: str, in_degree: int, max_clause: int = 20) -> str:
        max_clauses: int = int(
            factorial(in_degree) / (factorial(in_degree // 2)
                                    * factorial(in_degree - in_degree // 2))
        )
        return f'maxC("{n}", {min(max_clauses, max_clause)}).'


# ------------------------------------------------------------------------------
# Timeseries Observations description
# ------------------------------------------------------------------------------
class TimeSeries:
    # --------------------------------------------------------------------------
    # Template ASP
    # --------------------------------------------------------------------------
    @classmethod
    def constraint_bounds(cls, exp_id: str, r: str,
                          lb: float, ub: float) -> str:
        return f'param(transport, "{exp_id}", "{r}", "{lb}", "{ub}").'

    @classmethod
    def constraint_mutant(cls, exp_id: str, n: str, v: bool) -> str:
        return f'param(mutation, "{exp_id}", "{n}", "{v}").'

    @classmethod
    def constraint_datatype(cls, exp_id: str,
                            datatypes: set[Literal['kinetics', 'fluxomics',
                                                   'transcriptomics']]) -> str:
        datatypes_asp: list[str] = [
            f'datatype("{exp_id}", "{dt}").'
            for dt in sorted(datatypes)
        ]
        return '\n'.join(datatypes_asp)

    @classmethod
    def observation_value(cls, t: tuple[str, int], r: str, v: float,
                          threshold: float = 10**-8) -> str:
        v_bin: int = 1 if v > threshold else -1
        return f'obs(("{t[0]}", {t[1]}), "{r}", {v_bin}).'

    @classmethod
    def observation_growth(cls, t: tuple[str, int], v: float) -> str:
        return f'obj(("{t[0]}", {t[1]}), "{v}").'

    @classmethod
    def observation_transition(cls, t1: tuple[str, int], t2: tuple[str, int]) \
            -> str:
        return f'next(("{t1[0]}", {t1[1]}), ("{t2[0]}", {t2[1]})).'

    @classmethod
    def observation_exchange(cls, t: tuple[str, int], r: str,
                             lb: float, ub: float) -> str:
        return f'bound(("{t[0]}", {t[1]}), "{r}", "{lb}", "{ub}").'

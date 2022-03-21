#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Kerian Thuillier"
__email__ = "kerian.thuillier@ens-rennes.fr"

#IMPORT#########################################################################

from scipy.special import binom
from .Parameters import (
    ANNOTATION_COLUMN,
    DATA_TYPE_COLUMN,
    MAX_OBS_BUFFER,
    THRESHOLD,
    UpdateMode)
from .MetabolicNetworks import MetabolicNetwork

from .Parameters import BIOMASS

#MODEL##########################################################################

inferring_model = """
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRÉ-CALCUL SUR LES OBSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time(T1) :- next(T1, _).
time(T2) :- next(_, T2).

experiment(E) :- time((E, _)).

nb_obs(E, S) :- S = #count {T: time((E, T)) }, experiment(E).
1 { totalTimes(E, S..S+K) } 1 :- experiment(E), nb_obs(E, S), maxObsToAdd(K).

:- tlink(_, (E, TM)), totalTimes(E, S), S < TM.

ttime((E, TM)) :- totalTimes(E, S), TM = (1..S). 
tnext((E, TM1), (E, TM2)) :- ttime((E, TM1)), ttime((E, TM2)), TM1 + 1 = TM2.

:- next((E, TO1), (E, TO2)), tlink((E, TO1), (E, TM1)), tlink((E, TO2), (E, TM2)), TM2 <= TM1.
:- next((E, TO1), (E, TO2)), tlink((E, TO1), (E, TM1)), tlink((E, TO2), (E, TM2)), maxObsToAdd(K), K + 1 < TM2 - TM1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRÉPARATION DES DONNÉES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inp(X,R) :- reactant(X,R,_), not product(X,_,_).
out(X,R) :- product(X,R,_), not reactant(X,_,_).
r(r,A,R) :- reactant(A,R,_), product(A,_,_).
r(p,A,R) :- product(A,R,_), reactant(A,_,_).

varm(A) :- r(_,A,_).
varm(A) :- r(_,_,A).
varm(A) :- inp(A,_).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BONESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{clause(N,1..C,L,S): in(L,N,S), maxC(N,C), node(N), node(L)}.

:- clause(N,_,L,S), clause(N,_,L,-S).

1 {constant(N,(-1;1)) } 1 :- node(N), not clause(N,_,_,_).

constant(N) :- constant(N,_).

size(N,C,X) :- X = #count {L,S: clause(N,C,L,S)}; clause(N,C,_,_).

:- clause(N,C,_,_); not clause(N,C-1,_,_); C > 1.

:- size(N,C1,X1); size(N,C2,X2); X1 < X2; C1 > C2.

:- size(N,C1,X); size(N,C2,X); C1 > C2; 
    mindiff(N,C1,C2,L1) ; mindiff(N,C2,C1,L2) ; L1 < L2.

clausediff(N,C1,C2,L) :- 
    clause(N,C1,L,_);not clause(N,C2,L,_);clause(N,C2,_,_), C1 != C2.

mindiff(N,C1,C2,L) :- clausediff(N,C1,C2,L); 
    L <= L' : clausediff(N,C1,C2,L'), clause(N,C1,L',_), C1!=C2.

:- size(N,C1,X1); size(N,C2,X2); C1 != C2; X1 <= X2; 
    clause(N,C2,L,S) : clause(N,C1,L,S).

nbnode(NB) :- NB = #count{N: node(N)}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DÉFINITION DU MSS RESPECTANT LES OBSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ÉTAT MÉTABOLIQUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Définition des noeuds.
1 { v(T,A,(1;-1)) } 1 :- ttime(T), varm(A).

% Les noeuds doivent suivre les observations.
:- tlink(TO, TM), obs(TO,A,V), v(TM,A,-V).

% Un métabolite est produit/consommé par au moins une réaction.
:- ttime(T), r(S,A,_), v(T,A,1), v(T,R,-1): r(S,A,R).

% Une réactive active ses réactants et ses produits.
:- ttime(T), r(_,A,R), v(T,R,1), v(T,A,-1).

% Une réaction d'import doit avoir ses réactants dans l'environnement.
:- ttime(T), inp(X,R), v(T,X,-1), v(T,R,1).

1 { v(T,A,(1;-1)) } 1 :- varx(A), ttime(T).

% forward non-emss variables
:- varx(A), w(T,A,V), v(T,A,-V).

% if regulated is 0, mss cannot activate it
:- w(T,A,-1), v(T,A,1), node(A).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RÉSEAU BOOLÉEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% tune BoNesis encoding for using non-defined nodes
{clause(N,1..C,L,S): in(L,N,S), maxC(N,C), node(N)}.

read(T,A,V) :- tnext(T,_), not inp(A,_), v(T,A,V).
read(T,A,V) :- tnext(T,T2), inp(A,_), tlink(TO2, T2), obs(TO2,A,V).

%% eval
eval(T,A,C,-1)  :- update(T,A), clause(A,C,L,V), read(T,L,-V).
eval(T,A,C,1)   :- read(T,L,V): clause(A,C,L,V); update(T,A), clause(A,C,_,_).
eval(T,A,1)     :- eval(T,A,C,1), clause(A,C,_,_).
eval(T,A,-1)    :- eval(T,A,C,-1): clause(A,C,_,_); update(T,A), clause(A,C,_,_).
eval(T,A,V)     :- update(T,A), constant(A,V).

%% intermediate regulated state
mode(T1,reg) :- tnext(T1,_).

% copy inputs
w(T2,A,V) :- inp(A,_), tnext(_,T2), tlink(TO2, T2), obs(TO2,A,V).
% copy non-updated
w(T2,A,V) :- tnext(T1,T2), not inp(A,_), not update(T1,A), v(T1,A,V).
% apply update
w(T2,A,V) :- tnext(T1,T2), update(T1,A), eval(T1,A,V).

%% variables not in emss
varx(A) :- node(A), not varm(A).

% TODO: 
%   coder la vraie règle (peut passer à 0 seulement si un flux le consomme)?
%   ou autoriser tout changement (revenir au concept de 'control' nodes 
%   dont on ne cherche pas à justifier la valeur)
constant(A,-1) :- inp(A,_).

% no constant
:- constant(A), not inp(A,_).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DÉBOGAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propre a nos donnees issues de FlexFlux (a (X, 0) aucun flux n'est calcule)
v(TM1,R,-1) :- r(_,_,R), next(TO1,_), TO1=(_,0), tlink(TO1,TM1).

% For monotonic constraints
w(T,A,1) :- tnext(_, T), not node(A), r(_,_,A).
w_obs(TO, A, 1) :- next(_, TO), not node(A), r(_, _, A).
w_obs(TO, A, V) :- tlink(TO, TM), w(TM, A, V).
:- tlink(TO, TM), node(A), w_obs(TO, A, V), w(TM, A, -V).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS/OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
:- tnext(T1, T2), inp(X, R), v(T1, X, 1), v(T2, X,-1), v(T1, R,-1).
:- tnext(T1, T2), out(X, R), v(T1, X,-1), v(T2, X, 1), v(T1, R,-1).
"""


#INPUT#DATA#####################################################################

def asp_metabolic_network(metabolic_network: MetabolicNetwork) -> str:
    mn_asp = []
    stoichiometric_coeff = metabolic_network.get_stoichiometric_coefficients()
    for (r, m), v in stoichiometric_coeff.items():
        if v > 0:
            mn_asp.append(f'product(\"{m}\",\"{r}\",\"{v}\").')
        else:
            mn_asp.append(f'reactant(\"{m}\",\"{r}\",\"{v}\").')
    mn_asp = '\n'.join(mn_asp)
    return (mn_asp)


def asp_fixed_nodes(nodes: set):
    pkn_asp = []
    for n in nodes:
        pkn_asp.append(f'node(\"{n}\").')
    return '\n'.join(pkn_asp)


def asp_pkn(pkn, max_clause_size: int = -1) -> str:
    pkn_asp = []

    for n in pkn.nodes:
        pkn_asp.append(f'{{node(\"{n}\")}}.')

    for (f, t, s) in pkn.edges:
        if s == 1:
            pkn_asp.append(f'in(\"{f}\", \"{t}\",  1).')
        elif s == -1:
            pkn_asp.append(f'in(\"{f}\", \"{t}\", -1).')
        else:
            pkn_asp.append(f'in(\"{f}\", \"{t}\", (-1;1)).')

    for (n, i) in pkn.in_degree(pkn.nodes()):
        maxC = int(binom(i, i//2))
        if max_clause_size != -1:
            maxC = min(maxC, max_clause_size)
        pkn_asp.append(f'maxC(\"{n}\", {maxC}).')

    pkn_asp = '\n'.join(pkn_asp)

    return pkn_asp


def asp_observations(observations, max_obs_buffer: int = MAX_OBS_BUFFER) -> str:
    observations_asp = []
    times = observations.index.values.tolist()
    for (exp, t) in times:
        nodes = observations.loc[(exp, t)]
        for (n, v) in nodes.items():
            if n in [ANNOTATION_COLUMN, DATA_TYPE_COLUMN, BIOMASS] or v is None:
                continue
            if v > THRESHOLD:
                observations_asp.append(f'obs(({exp}, {t}), \"{n}\", 1).')
            else:
                observations_asp.append(f'obs(({exp}, {t}), \"{n}\", -1).')

    obs2model_asp = []
    counter = None
    for i in range(1, len(times)):
        prev_time = times[i-1]
        curr_time = times[i]

        if prev_time[0] == curr_time[0]:
            observations_asp.append(f'next({prev_time}, {curr_time}).')

        if i-1 == 0:
            obs2model_asp.append(f'tlink({prev_time}, ({prev_time[0]}, 1)).')
            obs2model_asp.append(
                f'1 {{ tlink({curr_time}, ({curr_time[0]}, 2..2 + K)) }} 1 :- maxObsToAdd(K).')
            counter = 3
        elif prev_time[0] == curr_time[0]:
            obs2model_asp.append(
                f'1 {{ tlink({curr_time}, ({curr_time[0]}, {counter}..{counter} + K)) }} 1 :- maxObsToAdd(K).')
            counter += 1
        else:
            obs2model_asp.append(f'tlink({curr_time}, ({curr_time[0]}, 1)).')
            counter = 2

    observations_asp.append(f'maxObsToAdd({max_obs_buffer}).')

    observations_asp = '\n'.join(observations_asp)

    observations_asp += '\n' + '\n'.join(obs2model_asp)

    return observations_asp


#OPTIMISATION###################################################################


def optimisation_time_buffer() -> str:
    return '#minimize { V, E : totalTimes(E, V) }.'


#MODE###########################################################################


def update_mode(update_mode: UpdateMode) -> str:
    if update_mode == UpdateMode.SYNCHRONOUS:
        return 'update(T1,A) :- mode(T1,reg), node(A), not inp(A,_).'
    else:
        return '1{update(T1,A): node(A), not inp(A,_)} :- mode(T1,reg).'


#SHOW###########################################################################


def show_projection_network() -> str:
    return """
#show.
#show clause/4.
    """


def show_projection_component(component_nodes: set) -> str:
    component_nodes_str = [f'"{n}"' for n in component_nodes]
    component_node_asp = '(' + '; '.join(component_nodes_str) + ')'
    return f"""
#show.
#show clause(N,C,A,V): clause(N,C,A,V), N={component_node_asp}.
    """


def show_projection_node(node: str) -> str:
    node_str = f'"{node}"'
    return f"""
#show.
#show clause(N,C,A,V): clause(N,C,A,V), N={node_str}.
    """

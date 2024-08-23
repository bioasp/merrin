#program base.
% ==============================================================================
% Regulatory functions
% ==============================================================================

% ------------------------------------------------------------------------------
% Definitions
% ------------------------------------------------------------------------------

% Definition
{clause(N,1..C,L,S)} :- in(L,N,S), maxC(N,C), node(N).

% Clause have the smallest valid number possible
:- clause(N,C,_,_); not clause(N,C-1,_,_); C > 1.

% Regulatory functions are monotone
:- clause(N,_,L,S), clause(N,_,L,-S).

% No clause is a subset of another one
:- size(N,C,X), size(N,C',X'), C != C', X <= X',
clause(N,C',L,S): clause(N,C,L,S).
% [Equivalent ?]
:- clause(N,C,_,_), clause(N,C',_,_), C != C',
clause(N,C',L,S): clause(N,C,L,S).

% Nodes that do not have functions are constant
:- node(N), not clause(N,_,_,_).

% ------------------------------------------------------------------------------
% Orders
% ------------------------------------------------------------------------------

% Sort by length
size(N,C,X) :- clause(N,C,_,_), X=#count{L,S: clause(N,C,L,S)}.
:- size(N,C,X), size(N,C',X'), C < C', X > X'.

% Sort by lexicographic order
clausediff(N,C,C',L) :- clause(N,C,L,_), not clause(N,C',L,_), clause(N,C',_,_).
mindiff(N,C,C',L) :- clausediff(N,C,C',L), L <= L': clausediff(N,C,C',L').
:- size(N,C,X), size(N,C',X), C < C',
mindiff(N,C,C',L), mindiff(N,C',C,L'), L > L'.

% ==============================================================================
% Boolean Networks
% ==============================================================================

% ------------------------------------------------------------------------------
% Update mode
% ------------------------------------------------------------------------------

% Synchronous
update(T,A) :- out(T,A,_), not param(mutation,E,A,_), T=(E,_), show(A), node(A).

% ------------------------------------------------------------------------------
% Input Boolean state
% ------------------------------------------------------------------------------

% Elements that are not external metabolites are read as their state at T.
read(T,A,V) :- update(T,N), inp(T,A,V), in(A,N,_).

% ------------------------------------------------------------------------------
% Evaluation of DNF formula
% ------------------------------------------------------------------------------

% Clause evaluation
eval(T,A,C, 1) :- update(T,A), clause(A,C,_,_), read(T,L,V): clause(A,C,L,V).
eval(T,A,C,-1) :- update(T,A), clause(A,C,L,V), read(T,L,-V).

% Formula evaluation
eval(T,A, 1) :-  update(T,A), clause(A,C,_,_), eval(T,A,C,1).
eval(T,A,-1) :-  update(T,A), eval(T,A,C,-1): clause(A,C,_,_).
eval(T,A, 1)  :- update(T,A), not node(A).
eval(T,A, 1)  :- out(T,A,_), show(A), not update(T,A).

% ------------------------------------------------------------------------------
% Output Boolean state
% ------------------------------------------------------------------------------

:- eval(T,A,V), out(T,A,-V), show(A).

% ==============================================================================
% Reduced formula - Break the symmetry
% /!\ Warning, these constraints can be quite costly /!\
% ==============================================================================

% ------------------------------------------------------------------------------
% Prohibit unnecessary literals in clause
% ------------------------------------------------------------------------------
valid_lit(T,N,C,A) :- clause(N,C,A,V), read(T,A,-V),
                      read(T,A',V'): clause(N,C,A',V'), A' != A.
valid_lit(T,N,C,A) :- read(T,A,_), clause(N,C,A,_), A = A': clause(N,C,A',_).
:- not valid_lit(_,N,C,A), clause(N,C,A,_).

% ------------------------------------------------------------------------------
% Prohibit unnecessary clauses
% ------------------------------------------------------------------------------
valid_clause(T,N,C) :- clause(N,C,_,_), eval(T,N,C,1),
                       eval(T,N,C',-1): clause(N,C',_,_), C != C'.
valid_clause(T,N,C) :- eval(T,N,C,_), clause(N,C,_,_),
                       C = C': clause(N,C',_,_).
:- not valid_clause(_,N,C), clause(N,C,_,_).

% ------------------------------------------------------------------------------
% Prohibit unnecessary literals in rules
% ------------------------------------------------------------------------------
valid_rule(N,C,A) :- valid_lit(T,N,C,A), valid_clause(T,N,C).
:- not valid_rule(N,C,A), clause(N,C,A,_).

% ==============================================================================
% Debug
% ==============================================================================
param(mutation,-1,-1,-1).

% ============================================================================ %
%                  MERRIN: Metabolic Regulatory Rules INference                %
% ============================================================================ %

% ==============================================================================
% Mapping observation times with simulation times
% ==============================================================================

% ------------------------------------------------------------------------------
% Simulation time
% ------------------------------------------------------------------------------

% Number of timestep by experiment
nbObs(E,S) :- obs((E,_),_,_), S=#count{I: obs(T,_,_), T=(E,I)}.

% Time definition
time((E,1..S)) :- obs((E,_),_,_), nbObs(E,S).

% Added time
1 { gap(E,0..K) } 1 :- obs((E,_),_,_), maxGap(K).
time((E,S+A)) :- obs((E,_),_,_), nbObs(E,S), gap(E,K), A=0..K.

% Simulation time must be successive
:- time((E,ID)), ID > 1, not time((E,ID-1)).

% Successive simulation time
maxTs(E,S+K) :- nbObs(E,S), gap(E,K).

% Successive time definition
succ((E,T1),(E,T2)) :- maxTs(E,ID), T1=(1..(ID-1)), T2=T1+1,
time((E,T1)), time((E,T2)).

% ------------------------------------------------------------------------------
% Linking observation time with simulation time
% ------------------------------------------------------------------------------

% Mapping between observations and simulation times
1 { map(To,Ts): time(Ts), Ts=(E,_) } 1 :- obs(To,_,_), To=(E,ID).

% Order must be preserved
:- map(To,Ts), map(To',Ts'), To=(E,ID), To'=(E,ID'), next(To,To'), Ts > Ts'.

% Bijective mapping
:- map(To,Ts), map(To',Ts'), To = To', Ts != Ts'.
:- map(To,Ts), map(To',Ts), To != To'.

% The first observation is the first simulation time
:- next(To,_), not next(_,To), Ts=(E,_), map(To,Ts), Ts!=(E,1).

% ==============================================================================
% Regulatory functions
% ==============================================================================

% ------------------------------------------------------------------------------
% Definitions
% ------------------------------------------------------------------------------

% Definition
{clause(N,1..C,L,S): in(L,N,S), maxC(N,C), node(N), node(L)}.
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
1 {constant(N,(-1;1))} 1 :- node(N), not clause(N,_,_,_).
constant(N) :- constant(N,_).

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

% ------------------------------------------------------------------------------
% Special case: reversible reactions clauses
% ------------------------------------------------------------------------------

% Reversible reactions have the same regulatory rules
node(Rr) :- rev(Rf,Rr), node(Rf).
node(Rf) :- rev(Rf,Rr), node(Rr).
clause(Rr,C,L,V) :- rev(Rf,Rr), clause(Rf,C,L,V).
clause(Rf,C,L,V) :- rev(Rf,Rr), clause(Rr,C,L,V).

% ==============================================================================
% Boolean Networks
% ==============================================================================

% ------------------------------------------------------------------------------
% Update mode
% ------------------------------------------------------------------------------

% Synchronous
update(T,A) :- time(T), node(A), not ext(A), not param(mutation,E,A,_), T=(E,_).

% ------------------------------------------------------------------------------
% Input Boolean state
% ------------------------------------------------------------------------------

% Elements that are not external metabolites are read as their state at T.
read(T,A,V) :- not ext(A), not gene(A), succ(T,_), v(T,A,V), in(A,_,_).
read(T,A,V) :- not ext(A), not gene(A), succ(T,_), v(T,A,V), in(_,A,_).

% External metabolites are read as their state at T+1.
read(T,A,V) :- ext(A), succ(T,T'), v(T',A,V), in(A,_,_).
read(T,A,V) :- ext(A), succ(T,T'), v(T',A,V), in(_,A,_).

% Genes are read as their state at T+1.
read(T,A,V) :- gene(A), succ(T,T'), v(T',A,V), in(A,_,_).
read(T,A,V) :- gene(A), succ(T,T'), v(T',A,V), in(_,A,_).

% ------------------------------------------------------------------------------
% Evaluation of DNF formula
% ------------------------------------------------------------------------------

% Clause evaluation
eval(T,A,C, 1) :- update(T,A), clause(A,C,_,_), read(T,L,V): clause(A,C,L,V).
eval(T,A,C,-1) :- update(T,A), clause(A,C,L,V), read(T,L,-V).

% Formula evaluation
eval(T,A, 1) :- update(T,A), clause(A,C,_,_), eval(T,A,C,1).
eval(T,A,-1) :- update(T,A), eval(T,A,C,-1): clause(A,C,_,_).
eval(T,A,V)  :- update(T,A), constant(A,V).

% ------------------------------------------------------------------------------
% Output Boolean state
% ------------------------------------------------------------------------------

% --------------------------------------
% Fixed states
% --------------------------------------

% Value for mutated node
x(T',A,V) :- param(mutation,E,A,V), succ(_,T'), T'=(E,_).

% Value of external metabolites
x(T',A,V) :- ext(A), not r(_,_,A), succ(T,T'), read(T,A,V).

% Value for not updated nodes
x(T',A,V) :- succ(T,T'), not update(T,A), not ext(A), not param(mutation,E,A,_),
                x(T,A,V), T'=(E,_).

% --------------------------------------
% Updated States
% --------------------------------------

% Case Everythings: slow evaluation
x(T',A,V) :- not ext(A), not r(_,_,A), succ(T,T'), update(T,A), eval(T,A,V).

% Case Reaction: instant evaluation of its regulatory rule
x(T',A,V) :- not ext(A), r(_,_,R), succ(T,T'), update(T,A), eval(T,A,V).

% --------------------------------------
% Default States
% --------------------------------------

% Reactions that are not node are active by default
x(T',A,1) :- succ(_,T'), not node(A), r(_,_,A), not param(mutation,E,A,_),
    T'=(E,_).

% Elements that are not node are active by default
x(T',A,1) :- succ(_,T'), not node(A), in(_,A,_), not param(mutation,E,A,_),
    T'=(E,_).

% ------------------------------------------------------------------------------
% Constant formula
% ------------------------------------------------------------------------------

% External metabolite are disabled by default
constant(A,-1) :- ext(A).

% Constant function are prohibited
:- constant(A), not ext(A).

% ==============================================================================
% Boolean metabolic steady states
% ==============================================================================

% Internal structure of the metabolic network
r(r,A,R) :- reactant(A,R,_), not ext(A).
r(p,A,R) :- product(A,R,_), not ext(A).

% All the elements of the metabolic network
varm(A) :- r(_,A,_).
varm(R) :- r(_,_,R).
varm(A) :- ext(A).

% All the elements of the regulatory network
varx(A) :- in(_,A,_), not varm(A).

% ------------------------------------------------------------------------------
% Definitions
% ------------------------------------------------------------------------------

% Binary states associated with each element
1 { v(T,A,(-1;1)) } 1 :- time(T), varm(A).
1 { v(T,A,(-1;1)) } 1 :- time(T), varx(A).

% ------------------------------------------------------------------------------
% Metabolic network properties
% ------------------------------------------------------------------------------

% An active internal metabolite must produced/consummed
:- r(S,A,_), v(T,A,1), v(T,R,-1): r(S,A,R).

% All products and reactants of an active internal reaction are active
:- r(_,A,R), v(T,R,1), v(T,A,-1).

% Import reaction has all its reactants at T
:- ext(A), reactant(A,R,_), v(T,R,1), v(T,A,-1).

% Reversible reaction can only be activated in one direction
:- rev(Rf,Rr), Rf!=-1, Rr!=-1, v(T,Rf,1), v(T,Rr,1).

% ------------------------------------------------------------------------------
% Regulations
% ------------------------------------------------------------------------------

% State of regulator depends of the regulatory state
:- varx(A), x(T,A,V), v(T,A,-V).

% Reaction can be inhibited by the regulatory state
:- r(_,_,R), node(R), x(T,R,-1), v(T,R,1).

% ------------------------------------------------------------------------------
% Necessarily activated reactions
% ------------------------------------------------------------------------------
:- fixedFlux(R1), r(_,_,R2), v(T,R1,-1), v(T,R2,1).

% ==============================================================================
% Successive Boolean Metabolic Steady-States Dynamics
% ==============================================================================

% ------------------------------------------------------------------------------
% External metabolite dynamics
% ------------------------------------------------------------------------------

% New external metabolites should be produced at T-1
:- succ(_,T), succ(T,T'), ext(A), v(T,A,-1), v(T',A,1), v(T,R,-1): product(A,R,_).

% Removed external metabolites should be consummed at T-1
:- succ(T,T'), ext(A), v(T,A,1), v(T',A,-1), v(T,R,-1): reactant(A,R,_).

% ==============================================================================
% Linear Constraints - Checking Boolean metabolic steady state optimality
%                      and ensure that the observed lower bound is satisfied
% ==============================================================================

% ------------------------------------------------------------------------------
% Primal FBA model
% ------------------------------------------------------------------------------

lpBound(Ts,R,0,0) :- time(Ts),     x(Ts,R,-1), r(_,_,R).
lpBound(Ts,R,0,0) :- time(Ts), not x(Ts,R,-1), reactant(M,R,_), ext(M), v(Ts,M,-1).
lpBound(Ts,R,L,U) :- time(Ts), not x(Ts,R,-1), v(Ts,M,1): ext(M), reactant(M,R,_); r(_,_,R), bound(R,L,U), Ts=(E,_), not param(transport,E,R,_,_), not map(_,Ts).
lpBound(Ts,R,L,U) :- time(Ts), not x(Ts,R,-1), v(Ts,M,1): ext(M), reactant(M,R,_); Ts=(E,_), param(transport,E,R,L,U), not map(_,Ts).
lpBound(Ts,R,L,U) :- time(Ts), not x(Ts,R,-1), v(Ts,M,1): ext(M), reactant(M,R,_); r(_,_,R), bound(R,L,U), Ts=(E,_), not param(transport,E,R,_,_), map(To,Ts), not bound(To,R,_,_).
lpBound(Ts,R,L,U) :- time(Ts), not x(Ts,R,-1), v(Ts,M,1): ext(M), reactant(M,R,_); Ts=(E,_), param(transport,E,R,L,U), map(To,Ts), not bound(To,R,_,_).
lpBound(Ts,R,L,U) :- time(Ts), not x(Ts,R,-1), v(Ts,M,1): ext(M), reactant(M,R,_); map(To,Ts), bound(To,R,L,U).

% Variables
&dom(check(Ts)){L..U} = f(R) :- lpBound(Ts,R,L,U).

% Steady-state
&sum(check(Ts)){S * f(R): reactant(M,R,S);
S * f(R): product(M,R,S)} = 0 :- time(Ts), not ext(M), r(_,M,_).

% Inhibition due to missing input metabolite in substrate
&sum(check(Ts)){f(R)} = 0 :- time(Ts), ext(M), reactant(M,R,_), v(Ts,M,-1).

% Inhibition due to regulatory rules
&sum(check(Ts)){f(R)} = 0 :- time(Ts), r(_,_,R), x(Ts,R,-1).

% Forced metabolic flux for reactions of interest
&sum(check(Ts)){f(R)}  = 0 :- time(Ts), in(R,_,_), r(_,_,R), v(Ts,R,-1).
&sum(check(Ts)){f(R)} >= E :- time(Ts), in(R,_,_), r(_,_,R), v(Ts,R, 1),
epsilon(E).

% If a fixed reaction is set to zero then the optimum flux should be 0
&assert(check(Ts)){f(Obj)} <= E :- lpcheck_ignore(Ts), objective(Obj),
                                   epsilon(E).
&sum(check(Ts)){f(Obj)} >= LB * O :- time(Ts), map(To,Ts), objective(Obj),
                                        lb(LB), obj(To,O).

% ------------------------------------------------------------------------------
% Dual FBA model - to ensure optimality
% ------------------------------------------------------------------------------
lpcheck_ignore(Ts) :- time(Ts), fixedFlux(R), R != -1, v(Ts,R,-1).
lpcheck_dual(Ts) :- time(Ts), not map(_,Ts).
lpcheck_dual(Ts) :- time(Ts), map(To,Ts), not obj(To,_).

% Variables
&sum(check(Ts)){p(R)} <= 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts),
                              r(_,_,R).
&sum(check(Ts)){q(R)} >= 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts),
                              r(_,_,R), not freeBound(Ts,R).
&sum(check(Ts)){s(R)} >= 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts),
                              r(_,_,R), freeBound(Ts,R).
&sum(check(Ts)){t(R)} >= 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts),
                              r(_,_,R), freeBound(Ts,R).

% Objective
freeBound(T1,R) :- reactant(M,R,_), ext(M), succ(T1,T2), v(T1,M,1), v(T2,M,-1).

&sum(check(Ts)){
    -dual;
    L * p(R): r(_,_,R), lpBound(Ts,R,L,_);
    U * q(R): r(_,_,R), lpBound(Ts,R,_,U), not freeBound(Ts,R);
    t(R): r(_,_,R), freeBound(Ts,R)
} = 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts).

% Constraints
&sum(check(Ts)){
    S * r(M): not ext(M), reactant(M,R,S);
    S * r(M): not ext(M), product(M,R,S);
    p(R);
    q(R): not freeBound(Ts,R);
    s(R): freeBound(Ts,R)
} = 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts),
         r(_,_,R), not objective(R).
&sum(check(Ts)){
    S * r(M): not ext(M), reactant(M,R,S);
    S * r(M): not ext(M), product(M,R,S);
    p(R);
    q(R): not freeBound(Ts,R);
    s(R): freeBound(Ts,R)
} = 1 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts), objective(R).

&sum(check(Ts)){
    -U * s(R);
    t(R)
} <= 0 :- lpcheck_dual(Ts), succ(_,Ts), not lpcheck_ignore(Ts), lpBound(Ts,R,_,U),
          freeBound(Ts,R).

% ------------------------------------------------------------------------------
% Ensure optimality
% ------------------------------------------------------------------------------
&sum(check(Ts)){ f(Obj); -dual} <= E :- lpcheck_dual(Ts), succ(_,Ts),
                                        not lpcheck_ignore(Ts), objective(Obj),
                                        epsilon(E).
&sum(check(Ts)){-f(Obj);  dual} <= E :- lpcheck_dual(Ts), succ(_,Ts),
                                        not lpcheck_ignore(Ts), objective(Obj),
                                        epsilon(E).

% ==============================================================================
% Linear Constraints - Ensure that no metabolic steady state has an optimum
%                      growth value than the observation
% ==============================================================================
lpMatchObjective(To) :- next(_,To), obj(To,_).

% Variables

% Metabolic flux bounds
&dom(reg(To)){L..U} = f(R) :- lpMatchObjective(To), param(transport,E,R,L,U),
To=(E,_).
&dom(reg(To)){L..U} = f(R) :- lpMatchObjective(To), bound(To,R,L,U).
&dom(reg(To)){L..U} = f(R) :- lpMatchObjective(To), r(_,_,R), bound(R,L,U).

% Steady-state
&sum(reg(To)){S * f(R): reactant(M,R,S);
S * f(R): product(M,R,S)} = 0 :- next(_,To),
not ext(M), r(_,M,_).

% Inhibition due to missing input metabolite in substrate
&sum(reg(To)){f(R)} = 0 :- next(_,To), map(To,Ts), ext(M),
reactant(M,R,_), v(Ts,M,-1).

% Inhibition due to regulatory rules
&sum(reg(To)){f(R)} = 0 :- next(_,To), map(To,Ts), r(_,_,R),
x(Ts,R,-1).

% Ensure that the biomass optimum match the observation
&assert(reg(To)){f(Obj)} <= UB * O :- next(_,To), objective(Obj),
   ub(UB), obj(To,O).

% ==============================================================================
% Applying observations
% ==============================================================================

% ------------------------------------------------------------------------------
% Boolean metabolic steady state
% ------------------------------------------------------------------------------

% Observations are fixed
:- obs(To,A,V), map(To,Ts), v(Ts,A,-V).

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
% Objective function
% ==============================================================================

% Minimise the number of simulation time to add
#minimize{S, E: gap(E,S)}.

% ==============================================================================
% Debug
% ==============================================================================

% Remove warnings for undefined optional elements
gene(-1).
obj(-1,-1).
rev(-1,-1).
fixedFlux(-1).
bound(-1,-1,-1).
bound(-1,-1,-1,-1).
param(transport,-1,-1,-1,-1).
param(mutation,-1,-1,-1).
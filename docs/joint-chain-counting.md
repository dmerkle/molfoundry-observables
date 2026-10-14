# Exact joint counting on the unbounded chain

Joint counting uses finite combinatorial partitions, with no assumption of
stochastic independence. The connected-footprint lemma and the supported local
condition semantics are prerequisites.

## 1. Cross-component NoEdge

For injective maps into a simple host, the complement of `NoEdge(u,v)` is a
disjoint union: precisely one positive edge `u --b-- v` for some label `b` in
the host language's finite edge signature. Therefore

    count(P, NoEdge(u,v) AND c)
      = count(P,c) - sum_b count(P + edge_b(u,v),c).

Only cross-component guards need this expansion; internal guards are already
handled by the connected footprint counter. Each positive edge joins two
components. A support-type proof can discharge guards that can never see a
host edge. All remaining guards and full-host neighbour conditions are retained.

## 2. Equality-kernel partition identity

Let P have positive connected components P_1,...,P_k and no unresolved cross
guards. The product of their injective counts enumerates maps that are
injective on each component, but may collide between components. Each such
map has a unique kernel partition of the pattern vertices. A kernel block
contains at most one vertex from each original component and equal labels.

For each nontrivial compatible partition pi, maps with exactly kernel pi are
in bijection with **globally injective** maps of the quotient P/pi. Conditions
are pushed to the identified centres/endpoints. Conflicting edge labels or
self-loops give no maps. NoEdge whose endpoints coincide is true in a loopless
host. A block's intersected support types may prove it impossible before
enumeration; this filter is conservative.

Thus

    inj(P) = product_i inj(P_i) - sum_{nontrivial kernels pi} inj(P/pi).

Each partition is subtracted exactly once; this is not a sum over single pair
collisions, and no alternating Möbius coefficient belongs in this identity.
Identically serialized quotient patterns are grouped with the number of distinct kernels
as their coefficient. Those multiplicities must not be deduplicated away.
The serializer renumbers vertices but does not compute complete graph-isomorphism
classes; additional isomorphism grouping remains a possible optimization.

Quotients have fewer vertices. Cross-guard expansion either removes a guard or
reduces the number of positive components. These measures establish finite
recursion. Every base case is a connected footprint count. Exact polynomial
plus finite-point sequence arithmetic preserves short-chain exceptions.

## 3. Free-vertex acceleration

If q_l isolated vertices of label l carry no adjacency predicates or incident
NoEdges, every embedding of the remaining core uses exactly v_l vertices of
that label. If N_l(n) is the host's total number of vertices of label l, then

    inj(P)(n) = inj(core)(n) * product_l (N_l(n)-v_l)_(q_l).

The falling factorial counts injective choices among unused host vertices.
This accelerates high-order moment patterns without prescribing a target
observable basis. It applies even when the remaining core carries supported
conditions. It does not factor expectations of random variables.

## 4. Scope and budgets

The method supports the existing one-port decorated-path grammar with finite
labels and local neighbour-count/NoEdge predicates. It does not extend to
arbitrary nested conditions, arbitrary graph grammars, or dangling-aware
vertex-deletion jump observables. Kernel enumeration may grow combinatorially;
`joint_budget` limits recursive/partition work and yields UNKNOWN on exhaustion,
never an exact partial count. Completed subqueries may safely remain cached.

Tests check known all-n path formulae (including n=1000 and order20), full-host
injective counts on randomized disconnected/conditional patterns, short-chain
point exceptions, and explicit fail-closed budget behavior.

This is an implementation extension using standard combinatorial partition
identities. It is not advertised as a new general graph-counting theorem.

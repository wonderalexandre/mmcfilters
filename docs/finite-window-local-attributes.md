# Finite-window local attributes

This document defines the public structural contract implemented by
`FiniteWindowLocalAttributeComputer`. The computation applies to a
`MorphologicalTree` over a finite two-dimensional image domain. It is shared by
component trees and trees of shapes; the public algorithm does not select a
different localization path from the descriptive tree kind.

## Observation window

An `ObservationWindow` is an ordered collection

\f[
W=(\delta_0,\ldots,\delta_{m-1})\subset\mathbb Z^2.
\f]

Its order is semantic: offset `j` always corresponds to bit `j` of a
`BinaryVisibilityState`. The C++ implementation requires:

- between 1 and 32 offsets;
- the zero offset `(0, 0)` exactly once;
- no duplicate offsets.

`WindowOffset` uses the fields `rowOffset` and `columnOffset`. A positive row
offset moves down and a positive column offset moves right. A translated sample
outside the image domain has no anchored entry; it is not mapped to the root or
to another live node.

## Anchored entries

For anchor pixel \(p\), let \(P(p)\) be its smallest node. The anchor branch is
the chain from \(P(p)\) to the root. For a valid translated sample \(q\), the
anchored entry is the first node on that branch whose support contains \(q\).
With \(A=P(p)\) and \(S=P(q)\), it is

\f[
E^p(q)=
\begin{cases}
A,&S\subseteq A,\\
S,&A\subset S,\\
\operatorname{LCA}(A,S),&A\text{ and }S\text{ are incomparable}.
\end{cases}
\f]

The LCA case is part of the generic connected-subset-tree formulation. In the
component-tree bitquad setting studied by the local-attribute method, relevant
nearby smallest nodes are locally comparable. This does not imply that every
pair drawn from every component-tree window is comparable. A tree of shapes
may require the LCA case even for nearby samples.

`anchoredEntryMap` preserves observation-window coordinates and represents a
missing entry with `std::nullopt`. `orderedAnchoredEntries` groups equal entries
into `AnchoredEntryMask` records and orders them by increasing inclusion, from
the anchor's smallest node toward the root. This branch order is distinct from
bottom-up aggregation over the whole tree.

## Visibility state and local rule

For a node \(X\) on the anchor branch, coordinate \(j\) of the binary
visibility state is one precisely when the anchored entry of \f$\delta_j\f$ is
contained in \(X\). State bits only change from zero to one while moving toward
the root, and all offsets sharing an anchored entry change simultaneously.

A `LocalRule` evaluates one `BinaryVisibilityState` into a value in an additive
Abelian group. Its C++ model provides:

```cpp
using Value = ...;

Value additiveIdentity() const;
Value evaluateLocalRule(BinaryVisibilityState state) const;
void addAssign(Value& target, const Value& source) const;
void subtractAssign(Value& target, const Value& source) const;
```

C++ verifies the presence and types of these operations. The rule author is
responsible for their algebraic laws: identity, associative and commutative
addition, and additive inverses. Arbitrary transition or merge policies are not
part of the public theorem-facing API.

## Explicit computation stages

The public generic C++ pipeline separates three roles:

1. `computeEventDeltas(tree, anchorPixel, window, rule)` returns one
   `EventDelta<Value>` for each distinct anchored entry of that anchor. Its
   fields identify the `anchorPixel`, the `anchoredEntry`, and the signed
   local-rule difference `value`. The first difference is measured from the
   additive identity; later differences compare consecutive visibility states.
2. `computeLocalAttributeIncrements(tree, window, rule)` sums the event deltas
   from every anchor into one `LocalAttributeIncrement<Value>` per dense node
   slot. An increment may be signed and is not yet a node attribute.
3. `aggregateLocalAttributeIncrements(tree, increments, rule)` adds fully
   accumulated children into their parents and returns one
   `NodeAttribute<Value>` per dense node slot.

`FiniteWindowLocalAttributeComputer::compute` performs those same stages and
returns the final `NodeAttribute<Value>` records. The increment path evaluates
one anchor at a time and streams its event-delta values directly into the
dense node increments. It reuses fixed-capacity scratch storage bounded by the
32-coordinate observation-window contract; it does not allocate entry maps,
ordered-entry sequences, or event-delta vectors per pixel. The internal dense
buffer stores only additive values while the computation is running; node
identifiers are materialized only in public normative result records. Concrete
multi-anchor rules, including bitquad rules, accumulate every anchor position
directly into one dense buffer instead of constructing and combining one full
temporary buffer per position. The four canonical bitquad windows and the
canonical contour-side window are immutable internal objects constructed once
and reused by subsequent computations. The explicit
`computeEventDeltas` operation still materializes and returns the normative
records when a caller requests that stage. Concrete bitquad and contour-side
storage remains an attribute-computer detail.

## Canonical bitquad specialization

A bitquad state uses row-major spatial coordinates:

| Coordinate | Position | Code weight |
|---|---|---:|
| `z0` | top-left | 1 |
| `z1` | top-right | 2 |
| `z2` | bottom-left | 4 |
| `z3` | bottom-right | 8 |

The displayed word is `z3z2z1z0` and the code is
`z0 + 2*z1 + 4*z2 + 8*z3`. Serialized and exported histograms use this exact
meaning. Each nonempty framed cell is assigned to the lowest-index visible bit,
so the four anchor positions count it exactly once. Samples outside the image
domain are background; an `R` by `C` image therefore has `(R+1)*(C+1)` framed
cells.

The five nonempty families are:

- `Q1 = {1, 2, 4, 8}`;
- `Q2 = {3, 5, 10, 12}`;
- `QD = {6, 9}`;
- `Q3 = {7, 11, 13, 14}`;
- `Q4 = {15}`.

Code zero has no `BitquadFamily`. It is not a sixth family. The bitquad pipeline
keeps each mathematical role distinct:

1. `computeNonemptyBitquadStateHistogramIncrements` returns signed 15-bin
   `NonemptyBitquadStateHistogramIncrement` values for codes 1 through 15.
2. `aggregateNonemptyBitquadStateHistogramIncrements` returns aggregated
   15-bin `NonemptyBitquadStateHistogram` values.
3. `materializeEmptyBitquadCount` derives bin zero by
   `h0 = (R+1)*(C+1) - sum(h1,...,h15)` and returns full
   `BitquadStateHistogram` values.
4. `computeBitquadFamilyIncrements` and `aggregateBitquadFamilyIncrements`
   likewise separate signed `BitquadFamilyIncrement` values from aggregated
   `BitquadFamilyCounts` values. Both types contain exactly `q1`, `q2`, `qd`,
   `q3`, and `q4`.

State and family counting describe node supports and do not choose a
connectivity convention. `BitquadAttributeProjection::materializeBitquadAttributes`
is the later, explicit stage that converts family counts to scalar attributes
under a `BitquadConnectivityPolicy`.

## C++ example

```cpp
#include <mmcfilters/localAttributes/FiniteWindowLocalAttributeComputer.hpp>

using namespace mmcfilters::local_attributes;

struct VisibleSampleCountRule {
    using Value = int;

    Value additiveIdentity() const { return 0; }
    Value evaluateLocalRule(BinaryVisibilityState state) const {
        return std::popcount(state.bits());
    }
    void addAssign(Value& target, const Value& source) const { target += source; }
    void subtractAssign(Value& target, const Value& source) const { target -= source; }
};

ObservationWindow window{{0, 0}, {0, 1}};
auto nodeAttributes = FiniteWindowLocalAttributeComputer::compute(
    tree, window, VisibleSampleCountRule{});
int rootValue = nodeAttributes[tree.root()].value;
```

The generic finite-window layer is currently a C++ extension API. Python users
consume its built-in attribute results through `Attribute`; generic local-rule,
event-delta, increment, and node-attribute records are not exposed as Python
objects.

# Editing API

This guide describes safe local mutations, staged structural edits, altitude
updates, rollback, and the lifetime of derived state.

## Ownership boundary

`MorphologicalTree` owns topology:

- dense internal `NodeId` slots;
- parent/child relations;
- direct proper-part ownership;
- optional regular 2D geometry and adjacency semantics.

`WeightedMorphologicalTree<T>` owns a private topology and a dense altitude
buffer. Its `topology()` accessor returns `const MorphologicalTree&`, so callers
cannot bypass the weighted edit boundary.

## Safe local edits

Two semantically complete mutations are public on `MorphologicalTree` and
`WeightedMorphologicalTree<T>`:

- `pruneNode(NodeId nodeId)`;
- `mergeNodeIntoParent(NodeId nodeId)`.

They reject invalid IDs, dead nodes, and the root before mutation. A successful
call publishes a committed tree, advances the topology mutation version, and
invalidates derived state computed against the previous version.

## Staged structural edits

An edit that temporarily detaches nodes, creates intermediate nodes, moves
children or proper parts, or changes the root must use an editor session:

```cpp
auto editor = tree.edit();
const NodeId inserted = editor.createDetachedNode();

editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);

editor.commit();
```

`MorphologicalTree::edit()` is the public factory; callers cannot construct a
`TreeEditor` directly. While an editor is open, the topology may be incomplete
and committed-tree operations are unavailable.

`MorphologicalTree` cannot be moved or cloned during an edit session. Finish the
session with `commit()` or `rollback()` before moving or cloning the tree. The
same rule applies to `WeightedMorphologicalTree<T>`.

### Commit and repair

`validateAndCommit()` validates one connected rooted tree, valid ownership, and
non-empty support for every live node. It returns `TreeValidationResult`:

- on success, it publishes the edit and closes the session;
- on failure, it returns `ok == false`, keeps the session open, and lets the
  caller repair the staged topology.

`commit()` is the exception-based wrapper. It closes the session on success and
throws when validation fails.

Complete validation is linear in the internal node slots plus proper parts:

```text
O(numInternalNodeSlots + numTotalProperParts)
```

There is no public unchecked commit path.

### Rollback

Public editors provide a strong rollback guarantee:

```cpp
auto editor = tree.edit();
editor.reparent(child, newParent);

if (shouldAbort) {
    editor.rollback();
} else {
    editor.commit();
}
```

`rollback()` restores the topology, ownership, slot state, root, and mutation
version from before the session. Destroying an active editor performs the same
rollback, including during exception unwinding.

### Incremental validation

For edits composed of supported primitives, advanced C++ callers can request a
move-only proof and publish the exact staged revision:

```cpp
auto editor = tree.edit();
// staged mutations

auto proof = editor.proveIncremental();
editor.commit(std::move(proof));
```

The proof is bound to one editor, tree, and mutation version. It cannot be copied
or reused after another mutation. If a primitive has no incremental validator,
proof construction uses complete validation.

## Weighted edits

`WeightedMorphologicalTree<T>::edit()` returns `WeightedTreeEditor<T>`, which
updates topology and altitude as one staged state:

```cpp
auto editor = weighted.edit();
const NodeId inserted = editor.createDetachedNode(insertedAltitude);

editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);

editor.commit();
```

`validateAndCommit()` validates topology first and then the declared altitude
order. A failed validation leaves the session open so the caller can repair the
topology or altitudes.

For ordered hierarchies:

- `INCREASING_FROM_ROOT` requires `altitude(parent) < altitude(child)`;
- `DECREASING_FROM_ROOT` requires `altitude(parent) > altitude(child)`;
- `UNCONSTRAINED` imposes no global direction but still requires valid finite
  values and a correctly sized buffer.

## Altitude setters

Committed weighted trees expose checked setters:

- `setAltitude(nodeId, value)` validates the ID, value, parent arc, and child
  arcs in `O(degree(nodeId))`;
- `setAltitudeBuffer(buffer)` validates buffer shape, finite floating-point
  values, and all ordered parent/child arcs in `O(numInternalNodeSlots)`.

Staged algorithms use `WeightedTreeEditor<T>::setNodeAltitude()` and cannot
publish until the weighted commit succeeds. Unchecked altitude setters are not
public API.

## Derived-state lifetime

Topology mutations invalidate objects that cache or interpret node-indexed
state. Guarded objects include:

- `ContoursComputedIncrementally::IncrementalContours`;
- `ContourTraceComputation::IncrementalContourTraces`;
- `WeightedTreeView<T>`;
- `AttributeFilters`;
- `ExtinctionValues`;
- `UltimateAttributeOpening`.

These objects capture the topology mutation version and reject reads when it no
longer matches.

Attribute buffers returned by value are not versioned. After topology mutation,
recompute attributes, contours, extinction values, and filter or
`UltimateAttributeOpening` helpers.

## Python boundary

Python exposes safe queries and local mutations:

- `pruneNode`;
- `mergeNodeIntoParent`;
- `setAltitude`;
- `setAltitudeBuffer` and `altitude`;
- queries such as `getRoot`, `getAliveNodeIds`, `getChildren`, and
  `getProperParts`.

Python does not expose `TreeEditor`, `WeightedTreeEditor`, `edit()`, unchecked
setters, or a mutable topology handle.

## Related guides

- [Morphological trees](trees.md): topology, ownership, and altitude contracts.
- [Attributes](attributes.md): recomputing node buffers after edits.
- [Filters](filters.md): helper lifetime after topology mutation.
- [Higra interoperability](higra-interoperability.md): imported-domain
  invalidation and fresh export snapshots.

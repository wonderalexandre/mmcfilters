# Editing API Contract

This project separates read-only topology, safe local mutations, and staged
structural edits. The goal is to keep ordinary callers away from low-level
topology rewiring while still allowing advanced algorithms to perform
multi-step edits explicitly.

For the broader tree ownership, altitude, and `NodeId` model, see
[Morphological Trees](trees.md).

## Public topology model

`MorphologicalTree` owns only the topology:

- dense internal `NodeId` storage;
- parent/child links;
- direct proper-part ownership;
- image-domain metadata and optional adjacency metadata.

Weighted values are not part of `MorphologicalTree`. They are owned by
`WeightedMorphologicalTree`, which stores:

- a private `MorphologicalTree` topology;
- a private dense altitude buffer indexed by internal `NodeId`.

`WeightedMorphologicalTree<T>::topology()` returns `const MorphologicalTree&`.
Callers can inspect the topology but cannot mutate it through that accessor.

## Safe committed edits

These methods are intentionally public on both `MorphologicalTree` and
`WeightedMorphologicalTree`:

- `pruneNode(NodeId nodeId)`
- `mergeNodeIntoParent(NodeId nodeId)`

They are local, semantically complete operations. They reject invalid ids,
dead nodes, and the root before mutating. They do not run the full connected
tree validation because they do not intentionally leave the tree in a staged
or disconnected state. They are still committed topology edits: they advance
the tree mutation version and invalidate topology-derived state computed before
the call.

For weighted trees, these methods delegate the topology change to the owned
`MorphologicalTree` and preserve the external altitude buffer.

## Staged structural edits

Any topology edit that can temporarily detach nodes, create intermediate
nodes, move children, move proper parts, or change the root must go through an
editor session:

```cpp
auto editor = tree.edit();
const NodeId inserted = editor.createDetachedNode();

editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);

editor.commit();
```

`TreeEditor` is not constructible by callers. `MorphologicalTree::edit()` is
the public factory. This keeps the edit boundary visible at each call site.

Opening an editor marks the tree as being in an edit session. The library does
not snapshot or roll back the previous state. If validation fails, the edit
session remains open and the caller can continue repairing the topology through
the same editor.

`TreeEditor::validateAndCommit()` runs the full connected-rooted-tree validation
and returns a `TreeValidationResult`. On success it closes the edit session. On
failure it returns `ok == false`, keeps the session open, and leaves the
partially edited topology in place.

`TreeEditor::commit()` is the exception-based wrapper around
`validateAndCommit()`: it closes the session on success and throws on failure.
Use `validateAndCommit()` when the caller wants to branch on a failed edit
without exception control flow; use `commit()` when failure should abort the
current operation.

The validation cost is linear in the current internal node slots plus proper
parts: `O(numInternalNodeSlots + numTotalProperParts)`.

For internal hot paths that preserve invariants by construction,
`TreeEditor::commitUnchecked()` closes the edit session without connected-tree
validation. It should not be used at public API boundaries. When
`MMCFILTERS_ENABLE_ASSERTS` is enabled, `commitUnchecked()` asserts that the
tree is valid before closing the session.

Destroying an active editor does not roll back or close the session. The tree
remains in editing mode, rejects committed-tree APIs, and rejects a new editor.
This makes an abandoned edit visible instead of silently publishing a partially
edited topology.

## Weighted staged edits

`WeightedMorphologicalTree<T>::edit()` returns `WeightedTreeEditor<T>`.

`WeightedTreeEditor<T>` wraps a structural `TreeEditor` and adds altitude-buffer
updates for new nodes:

```cpp
auto editor = weighted.edit();
const NodeId inserted = editor.createDetachedNode(insertedAltitude);

editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);

editor.commit();
```

`WeightedTreeEditor<T>::validateAndCommit()` first validates the topology and then
validates altitude order through
`WeightedMorphologicalTree<T>::validateMonotoneAltitude()`. If either check fails,
the weighted edit session remains open so the caller can repair topology or
altitudes. `WeightedTreeEditor<T>::commit()` is the exception-based wrapper.

`WeightedTreeEditor<T>::commitUnchecked()` is reserved for internal hot paths and
skips both topology and altitude-order validation.

For max-trees, altitude must be non-decreasing from parent to child. For
min-trees, altitude must be non-increasing from parent to child. Trees of
shapes currently skip monotone altitude validation.

## Altitude setters

Public altitude setters preserve the same committed-tree invariant:

- `WeightedMorphologicalTree<T>::setAltitude(nodeId, value)` validates the node id,
  the value domain, and the local altitude order against the node's parent and
  direct children before publishing the value. This local monotonicity check is
  `O(degree(nodeId))`.
- `WeightedMorphologicalTree<T>::setAltitudeBuffer(buffer)` validates the dense
  buffer shape, finite floating-point values, and full max-tree/min-tree
  monotonicity before replacing the owned buffer. The full monotonicity check is
  `O(numInternalNodeSlots)`.

`setAltitudeUnchecked()` and `setAltitudeBufferUnchecked()` are C++ escape
hatches for code that has already established altitude order by construction.
They still validate the committed-edit boundary, live node or buffer shape, and
finite floating-point values, but they intentionally skip monotone-altitude
validation. They are not exposed in Python.

## Derived-state lifetime

Tree topology mutations invalidate objects that cache or interpret node-indexed
state against a specific topology. This includes both staged edits published by
`TreeEditor` and safe committed edits such as `pruneNode()` and
`mergeNodeIntoParent()`. `MorphologicalTree` exposes a monotonic mutation
version for these guards. Objects that keep a reference to a tree capture that
version at construction and reject public reads after the tree changes.

The guarded objects include:

- `ContoursComputedIncrementally::IncrementalContours`
- `AttributeFilters`
- `ExtinctionValues`
- `UltimateAttributeOpening`

Plain attribute buffers returned by value are not versioned. After a topology
mutation, callers should recompute attributes, contours, extinction values, and
filter/UAO helper objects from the mutated tree.

## Python boundary

The Python API exposes the stable query surface and the safe mutators:

- `pruneNode`
- `mergeNodeIntoParent`
- `setAltitude`
- `setAltitudeBuffer` / `altitude`
- topology queries such as `getRoot`, `getAliveNodeIds`, `getChildren`, and
  `getProperParts`

Python does not expose `TreeEditor`, `WeightedTreeEditor`, `edit()`, or a
mutable topology handle from `WeightedMorphologicalTree`. It also does not
expose the C++ `Unchecked` altitude setters. If Python editor bindings are
added later, they should preserve the same commit boundary used by the C++ API.

## Regression guard

The contract is guarded by:

- `unit_edit_api_contracts`, which uses compile-time checks to prevent editor
  construction and low-level mutator exposure from becoming public;
- `unit_python_nodeid_api`, which ensures the Python API does not expose a
  mutable weighted topology handle or low-level structural mutators.

See `examples/editing_api_example.cpp` for a minimal compiled example.

## Related Guides

- [Morphological Trees](trees.md): owner/view boundary, altitude contracts, and
  mutation-version semantics.
- [Attributes](attributes.md): recomputing node-indexed buffers after topology
  edits.
- [Attribute Filters, Extinction Values, And UAO](filters.md): filter helper
  lifetime after tree mutation.
- [Higra Interoperability](higra-interoperability.md): imported Higra id-space
  invalidation and fresh export snapshots.

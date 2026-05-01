# Editing API Contract

This project separates read-only topology, safe local mutations, and staged
structural edits. The goal is to keep ordinary callers away from low-level
topology rewiring while still allowing advanced algorithms to perform
multi-step edits explicitly.

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

`WeightedMorphologicalTree::topology()` returns `const MorphologicalTree&`.
Callers can inspect the topology but cannot mutate it through that accessor.

## Safe public mutators

These methods are intentionally public on both `MorphologicalTree` and
`WeightedMorphologicalTree`:

- `pruneNode(NodeId nodeId)`
- `mergeNodeIntoParent(NodeId nodeId)`

They are local, semantically complete operations. They reject invalid ids,
dead nodes, and the root before mutating. They do not run the full connected
tree validation because they do not intentionally leave the tree in a staged
or disconnected state.

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

`TreeEditor::commit()` runs the full connected-rooted-tree validation. Its
dominant cost is linear in the current internal node slots plus proper parts:
`O(numInternalNodeSlots + numTotalProperParts)`.

## Weighted staged edits

`WeightedMorphologicalTree::edit()` returns `WeightedTreeEditor`.

`WeightedTreeEditor` wraps a structural `TreeEditor` and adds altitude-buffer
updates for new nodes:

```cpp
auto editor = weighted.edit();
const NodeId inserted = editor.createDetachedNode(insertedAltitude);

editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);

editor.commit();
```

`WeightedTreeEditor::commit()` first validates the topology through
`TreeEditor::commit()` and then validates the altitude order through
`WeightedMorphologicalTree::validateMonotoneAltitude()`.

For max-trees, altitude must be non-decreasing from parent to child. For
min-trees, altitude must be non-increasing from parent to child. Trees of
shapes currently skip monotone altitude validation.

## Python boundary

The Python API exposes the stable query surface and the safe mutators:

- `pruneNode`
- `mergeNodeIntoParent`
- topology queries such as `getRoot`, `getAliveNodeIds`, `getChildren`, and
  `getProperParts`

Python does not expose `TreeEditor`, `WeightedTreeEditor`, `edit()`, or a
mutable topology handle from `WeightedMorphologicalTree`. If Python editor
bindings are added later, they should preserve the same commit boundary used by
the C++ API.

## Regression guard

The contract is guarded by:

- `unit_edit_api_contracts`, which uses compile-time checks to prevent editor
  construction and low-level mutator exposure from becoming public;
- `unit_python_nodeid_api`, which ensures the Python API does not expose a
  mutable weighted topology handle or low-level structural mutators.

See `examples/editing_api_example.cpp` for a minimal compiled example.

# Attribute computer architecture

This guide documents the attribute-computer architecture and the extension path
for adding or changing attributes. For public usage, see [Attributes](attributes.md);
for the attribute table, see [Attribute catalog](attribute-catalog.md).

The public API calls each scalar an attribute. In this guide, descriptor refers
only to the mathematical quantity represented by a public scalar attribute.

In this subsystem, incremental means bottom-up and finite-window-oriented
computation over the current tree. It does not mean that every public attribute
buffer stays live after arbitrary topology edits.

## Public boundary

Ordinary application code should include `mmcfilters/attributes/Attributes.hpp`
and call `AttributeComputation`. Concrete computers are advanced extension
components, not alternate public orchestration APIs.

The C++ library is header-only, so installed packages include some `detail`
headers as transitive implementation dependencies. Those headers are shipped so
public headers compile downstream; they are not compatibility-contract headers.
The explicit installation manifest is
`cmake/mmcfiltersPublicHeaders.cmake`; adding a repository header does not
implicitly publish or install it.

## Computer contract

An attribute computer owns one coherent attribute family. Every computer must
provide:

- `inline static constexpr familyName` for diagnostics;
- `inline static constexpr family` for scheduler grouping;
- `inline static constexpr domain` for execution routing;
- `inline static constexpr producedAttributes` as the canonical attribute list;
- `static compute(context)` for internal-node rows;
- `static computeUnitRows(unitContext)` for compact exported-Higra unit rows.

Computers are stateless static kernels. The produced-attribute list has a single
source of truth: the computer class. `runtimeProducedAttributes<Computer>()`
materializes it only for call sites that need runtime storage.

Unit-row support is mandatory. If an attribute has a degenerate one-pixel
meaning, the computer defines that value explicitly. Otherwise it still defines
the exported unit-row convention.

## Registry and metadata

`AttributeRegistry.hpp` stores public attribute metadata:

- public name;
- description;
- group membership;
- altitude requirement;
- topology/support eligibility (`topologyOnly` in the registry).

Group membership is metadata. Public requests may mix scalar attributes and
groups; the pipeline expands groups, deduplicates scalars, and returns only the
requested public attributes.

`AttributeComputerRegistry.hpp` defines the computer protocol and
`RegisteredAttributeComputers`. Produced attributes are declared only by
`Computer::producedAttributes`, and scheduler grouping is declared only by
`Computer::family`.

## Execution model

At a high level, a request follows this path:

```text
request -> expand groups -> validate support -> materialize dependencies
        -> compute buffers -> assemble requested result -> project if needed
```

`AttributeFamilyScheduler` adds hidden dependencies, groups attributes by
family, and preserves dependency order. The central executors are:

- `executeAttributeComputationPlan(...)` for altitude-aware requests;
- `executeTopologyAttributeComputationPlan(...)` for topology/support requests.

The internal orchestration path is `detail::AttributePipeline`;
topology/support families are delegated to `TopologyAttributeBackend`. New code
should extend this path instead of adding another top-level execution pipeline.

Dependencies are ordinary attribute results consumed by another computer. They
are passed as `DependencySourceT<Real>`, a non-owning pair of `AttributeNames`
and `const Real*`. Dependency buffers are reusable only when they contain the
requested attributes and use `NodeIdSpace::MorphologicalTree`.

Several computers use bottom-up accumulation: preprocess the current node, merge
children into the parent, then finalize the current node. Increment-augmented public
calls compute the base attribute first, then materialize ancestor/descendant
sample offsets from a typed positive altitude step, sampling radius,
representative-descendant policy, and missing-sample policy.

## Contexts and concepts

The context types in `AttributeKernelSupport.hpp` are the adapter boundary:

- topology/support node rows use `AttributeComputeContext<Real>`;
- altitude-aware node rows use `AltitudeAttributeComputeContext<Real, T>`;
- topology/support unit rows use `UnitAttributeComputeContext<Real>`;
- altitude-aware unit rows use `AltitudeUnitAttributeComputeContext<Real, T>`.

`TopologyAttributeComputer` and `AltitudeAttributeComputer` enforce the standard
computer protocol. A new family should not add public family-specific method
names. Private helpers and `detail` kernels may keep narrower signatures when
that makes implementation or testing clearer.

The generic finite-window extension API exposes role-typed `EventDelta`,
`LocalAttributeIncrement`, and `NodeAttribute` records so its mathematical
pipeline can be inspected and tested. Attribute-family storage such as
`detail::BitquadFamilyCounts` and `detail::ContourSideCounts` remains an
implementation detail and is not part of the public attribute-computer
contract.

Bitquad-family counting is hierarchy-independent. Connectivity-dependent scalar
formulas are a later materialization step and receive an explicit
`detail::BitquadConnectivityPolicy`. For a tree of shapes with unequal
complementary connectivity, the policy consumes `ShapePolarity::Lower` or
`ShapePolarity::Upper` derived from exact node-versus-parent altitudes. The root
has no polarity and is handled by the policy's separate root entry.

## Numeric policy

Computers use `AttributeNumericPolicy.hpp` for degenerate divisions, square
roots, non-negative clamping, finite fallbacks, and ratio bounds. Attribute
buffers should not expose accidental `NaN` or infinite values for ordinary
finite inputs.

## Adding or changing attributes

Start by deciding whether the attribute belongs to an existing family or
requires a new family. Prefer an existing family when traversal, dependencies,
or intermediate state are shared.

Common metadata steps:

1. Add the scalar enum in `AttributeTypes.hpp` and one matching row in
   `AttributeRegistry.hpp`.
2. Classify the attribute as topology/support, altitude-aware,
   adjacency-dependent, or tree-kind specific.
3. Add it to a group only when the group semantics still hold.

For an attribute in an existing family:

1. Add it to the family's `producedAttributes`.
2. Extend request selection and `compute(context)`.
3. Use `DependencyResolver<Real>` for semantic dependencies and
   `AttributeNumericPolicy.hpp` for finite fallbacks.
4. Add attribute-level dependencies in `AttributeFamilyScheduler.hpp` only when
   another materialized attribute is consumed.
5. Extend `computeUnitRows(unitContext)`.
6. Add focused value tests, plus plumbing tests when registry, dependencies,
   projection, or public layout changes.

For a new family:

1. Add a new `AttributeComputerFamily` value.
2. Create a computer under `mmcfilters/attributes/computers/`.
3. Declare `familyName`, `family`, `domain`, and `producedAttributes`.
4. Implement `compute(context)` and `computeUnitRows(unitContext)`.
5. Register the computer in `RegisteredAttributeComputers`.
6. Register execution in `AttributePipeline.hpp` or `TopologyAttributeBackend.hpp`.
7. Add contract/plumbing tests and focused value tests.

When the public surface or attribute semantics change, update Python bindings
as needed and keep [Attributes](attributes.md) and the
[Attribute catalog](attribute-catalog.md) synchronized with the registry.

## Validation

Useful checks while changing this subsystem are:

```bash
cmake --build build --target \
  unit_public_attribute_api \
  unit_attribute_plumbing \
  unit_attribute_unit_values \
  unit_attributes_on_morphological_tree \
  unit_finite_window_local_attribute_computations \
  unit_maxdist_support

ctest --test-dir build --output-on-failure -R \
  "unit_(public_attribute_api|attribute_plumbing|attribute_unit_values|attributes_on_morphological_tree|finite_window_local_attribute_computations|maxdist_support|installed_consumer)"
```

Run Python tests when bindings or the Python facade change.

## Non-goals

- Do not introduce virtual attribute-computer classes.
- Do not move away from dense node ID buffers.
- Do not add runtime polymorphism to hot attribute kernels.
- Do not expose attribute-family-specific finite-window storage as public
  computer API; use the generic role-typed pipeline for extension work.

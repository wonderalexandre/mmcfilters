# Attribute Computer Architecture

This note records the current attribute-computer architecture. It is not a
transition plan.

## Contract

An attribute computer owns one coherent descriptor family. Every computer must
provide:

- `inline static constexpr familyName` for diagnostics;
- `inline static constexpr domain` for execution routing;
- `inline static constexpr producedAttributes` as the canonical descriptor list;
- `static compute(context)` for internal-node rows;
- `static computeUnitRows(unitContext)` for compact exported-Higra unit rows.

The produced-attribute list has a single source of truth: the computer class.
`runtimeProducedAttributes<Computer>()` materializes it as a
`std::vector<Attribute>` only for call sites that need runtime storage.

Unit-row support is mandatory. If a descriptor has a degenerate one-pixel
meaning, the computer defines that value explicitly. If the descriptor cannot
mathematically use internal-node state in the unit domain, the computer still
defines the exported unit-row convention.

## Protocol

`AttributeComputerProtocol.hpp` defines the common protocol:

- `AttributeComputerDomain`;
- `AttributeComputer`;
- `TopologyAttributeComputer`;
- `AltitudeAttributeComputer`;
- helpers such as `producesAttribute<Computer>(...)` and
  `runtimeProducedAttributes<Computer>()`;
- `RegisteredAttributeComputers`.

`Computer::domain` determines the compute context: topology families receive
topology contexts, while altitude families receive altitude-aware contexts.

Produced descriptors are declared only by `Computer::producedAttributes`.
Precise descriptor-level dependencies live in `AttributeFamilyScheduler.hpp`,
where the scheduler builds the recursive dependency closure for each public
request.

## Contexts

The context types in `AttributeKernelSupport.hpp` are the adapter boundary:

- topology/support node rows use `AttributeComputeContext<Real>`;
- altitude-aware node rows use `AltitudeAttributeComputeContext<Real, T>`;
- topology/support unit rows use `UnitAttributeComputeContext<Real>`;
- altitude-aware unit rows use `AltitudeUnitAttributeComputeContext<Real, T>`.

The contexts are borrowed views over topology, altitude spans, output buffers,
attribute layouts, requested subsets, and dependency sources. They do not own
storage.

## Concepts

`TopologyAttributeComputer` and `AltitudeAttributeComputer` enforce the standard
computer protocol. A new family should not add public family-specific method
names. Private helpers and `detail` kernels may keep narrow span-based
signatures when that makes the implementation clearer or easier to test.

## Execution

`AttributeFamilyScheduler` expands public requests, adds hidden dependencies,
groups descriptors by family, and preserves dependency order. The central
executors are:

- `executeAttributeComputationPlan(...)` for altitude-aware requests;
- `executeTopologyAttributeComputationPlan(...)` for topology/support requests.

The scheduler also owns the internal mapping between registered computers and
`AttributeFamily` ids. `familyForAttribute(...)` and unit-row projection derive
their dispatch from `RegisteredAttributeComputers`, so adding a family should
update the registered-computer list and the computer-to-family mapping together.

Public result layouts contain only requested attributes. Hidden dependencies
remain internal scratch data unless they were explicitly requested.

## Numeric Policy

Computers use `AttributeNumericPolicy.hpp` for degenerate divisions, square
roots, non-negative clamping, finite fallbacks, and ratio bounds. Attribute
buffers should not expose accidental `NaN` or infinite values for ordinary
finite inputs.

## Non-Goals

- Do not introduce virtual attribute-computer classes.
- Do not move away from dense node-id buffers.
- Do not add runtime polymorphism to hot attribute kernels.
- Do not expose local-event bucket storage as public computer API.

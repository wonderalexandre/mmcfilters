# Validation boundaries

MorphologicalAttributeFilters selects caller-contract validation at compile
time without changing the public C++ API. Normal builds use `CHECKED`;
controlled scientific benchmarks may use `UNCHECKED` when the experiment has
already established every input domain.

## Classification

Every condition belongs to one of the following categories:

| Category | Examples | Policy |
| --- | --- | --- |
| Caller contract | valid node or pixel id, non-null pointer, buffer shape, finite threshold, positive delta, supported padding, valid call sequence | `MMCFILTERS_CONTRACT_REQUIRE` or `MMCFILTERS_CONTRACT_CHECKED_ONLY` |
| Scientific capability | globally monotone altitude order required by MSER/extinction, adjacency information required by a projection | Always enforced when the algorithm cannot produce the defined scientific quantity without it |
| Structural invariant | coherent parent/child topology, valid support ownership, monotone owned altitude after mutation, projection conservation, representable output | Always enforced |
| Established kernel precondition | live-node domain, valid parent/child links, dense buffer indexing, compatible image domain already established by the caller | No runtime validation inside the kernel |

`UNCHECKED` removes only caller-contract diagnostics. Violating one of those
preconditions may cause undefined behavior. Scientific capability checks and
structural invariants do not become optional merely because a benchmark build
was requested.

## Control-flow rule

A public entry point validates its complete caller-owned input once and then
enters a validation-free execution core. It must not call another public
operation from the per-node or per-pixel loop.

Validation-free cores live in `detail::kernel`. They use established spans and
the narrow `CommittedTreeAccess`/`CommittedImageAccess` facades for primitive
tree and image access. Helpers outside `detail::kernel` may also use those
facades after their own dominating boundary has established the same domain.

For example, an MSER computation validates the altitude buffer, delta and
hierarchy capability once. Its ancestor search then reads altitude slots and
parent links directly for every live node. Likewise, a reconstruction filter
validates its criterion and output image once before traversing committed child
and proper-part ranges.

## Redundant checks

A defensive check is removed rather than wrapped when its truth follows from a
dominating boundary or from an always-maintained structural invariant. Typical
examples are rechecking that every node returned by the live-node iterator is
alive, revalidating an altitude span inside every ancestor search, and checking
the same topology mutation version in both an adapter and the implementation it
immediately calls.

The validation-boundary audit rejects contract macros, exceptions and defensive
public accessors inside `detail::kernel`. Tests run in both compile-time modes:
invalid-input expectations are enabled only in `CHECKED`, while scientific
result checks and invariant tests run in both modes.

## Benchmark rule

All algorithms in one comparison must be built with the same contract mode,
compiler options and optimization level. A benchmark should report the selected
mode and verify an identical scientific checksum between `CHECKED` and
`UNCHECKED` builds before timing results are interpreted.

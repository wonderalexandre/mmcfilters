# Scientific API validation-contract campaign

## Reproducibility

- Experiment Git commit: `84b5d415cebd6d12d8249a2467359d910b4a9164`
- Analysis Git commit: `84b5d415cebd6d12d8249a2467359d910b4a9164`; dirty worktree: `False`
- CPU: Apple M4
- Operating system: macOS-26.5.2-arm64-arm-64bit-Mach-O
- Compiler: Apple LLVM 21.0.0 (clang-2100.1.1.101)
- Profile: `publication`
- Repetitions per scenario/process: 15
- Process runs per contract mode: 5

Every experiment passed deterministic repetition checks and CHECKED/UNCHECKED checksum equality for all scenarios and structural outcomes.

## Workloads

| Workload | Domain | Scenarios | CHECKED sum (ms) | UNCHECKED sum (ms) | Difference |
| --- | --- | --- | --- | --- | --- |
| publication_structured | 512x512 | 141 | 13402.183 | 12342.681 | -7.91% |
| publication_real_lena | 256x256 | 145 | 1557.814 | 1526.714 | -2.00% |
| publication_real_brain | 322x506 | 145 | 2524.980 | 2510.670 | -0.57% |
| publication_real_wrist | 289x373 | 145 | 2370.015 | 2320.915 | -2.07% |

Across 576 workload-scenario observations, UNCHECKED was faster in 307. The median observation-level difference was -0.37%.

The descriptive sum of all independent scenario medians was 19854.992 ms in CHECKED and 18700.980 ms in UNCHECKED (-5.81%). It is not an application runtime.

## Raw-sample inference

Across the campaign, process-cluster bootstrap intervals support UNCHECKED as faster in 107 workload-scenario observations and slower in 34; 435 remain inconclusive.

| Workload | Scenario | Difference | Bootstrap 95% CI | Cliff's delta |
| --- | --- | --- | --- | --- |
| publication_structured | interoperability.established_input.normalized_altitude_saliency | -30.21% | [-33.68%, -21.93%] | -1.00 |
| publication_real_wrist | interoperability.established_input.normalized_altitude_saliency | -24.43% | [-27.95%, -13.01%] | -1.00 |
| publication_real_lena | interoperability.established_input.normalized_altitude_saliency | -23.72% | [-25.93%, -22.50%] | -1.00 |
| publication_real_brain | interoperability.established_input.normalized_altitude_saliency | -23.25% | [-23.76%, -14.72%] | -1.00 |
| publication_structured | editing.established_input.safe_merge_leaf | -23.25% | [-47.03%, -9.15%] | -1.00 |
| publication_real_lena | casf.established_input.incremental_area_medium_step_after_light | -21.57% | [-24.82%, -17.98%] | -1.00 |
| publication_structured | interoperability.established_input.topological_saliency | -21.40% | [-26.19%, -16.14%] | -1.00 |
| publication_structured | casf.established_input.incremental_sequence_int32_area | -21.38% | [-23.14%, -17.91%] | -1.00 |
| publication_structured | casf.established_input.incremental_area_heavy_step_after_medium | -21.29% | [-23.77%, -14.74%] | -1.00 |
| publication_real_wrist | interoperability.established_input.reconstruct_max_tree | -21.16% | [-23.21%, -19.83%] | -1.00 |
| publication_real_lena | interoperability.established_input.reconstruct_max_tree | -21.07% | [-22.07%, -19.73%] | -1.00 |
| publication_structured | interoperability.established_input.reconstruct_max_tree | -20.87% | [-21.96%, -16.99%] | -1.00 |

The campaign table preserves each workload-scenario interval separately; it does not pool inner repetitions or manufacture a larger independent sample.

## Suites and timing scopes

| Suite | Scope | Observations | CHECKED sum (ms) | UNCHECKED sum (ms) | Difference |
| --- | --- | --- | --- | --- | --- |
| attribute_bundles | established_input | 36 | 927.914 | 941.503 | +1.46% |
| attribute_groups | established_input | 48 | 3552.880 | 3649.666 | +2.72% |
| attributes | established_input | 264 | 1347.559 | 1374.594 | +2.01% |
| casf | end_to_end | 16 | 823.196 | 770.966 | -6.34% |
| casf | established_input | 28 | 4043.993 | 3278.721 | -18.92% |
| construction | end_to_end | 48 | 6525.391 | 6182.293 | -5.26% |
| editing | established_input | 16 | 5.523 | 5.239 | -5.13% |
| filters | established_input | 64 | 559.054 | 496.958 | -11.11% |
| interoperability | established_input | 36 | 771.278 | 728.692 | -5.52% |
| pipelines | end_to_end | 20 | 1298.205 | 1272.348 | -1.99% |

`end_to_end` and `established_input` contain different scientific operations. Their totals are descriptive and are not paired estimates of validation overhead.

## Paired end-to-end and established-input CASF

These rows have identical final checksums and threshold sequences. Their difference isolates CASF construction/setup from the established incremental sequence.

| Workload | Mode | End-to-end (ms) | Established (ms) | Setup (ms) | Setup share |
| --- | --- | --- | --- | --- | --- |
| publication_structured | CHECKED | 296.870 | 243.813 | 53.057 | 17.87% |
| publication_structured | UNCHECKED | 248.563 | 203.133 | 45.430 | 18.28% |
| publication_real_lena | CHECKED | 18.535 | 7.600 | 10.934 | 58.99% |
| publication_real_lena | UNCHECKED | 17.297 | 6.329 | 10.969 | 63.41% |
| publication_real_brain | CHECKED | 31.199 | 11.342 | 19.857 | 63.65% |
| publication_real_brain | UNCHECKED | 29.790 | 9.432 | 20.358 | 68.34% |
| publication_real_wrist | CHECKED | 30.024 | 13.200 | 16.824 | 56.03% |
| publication_real_wrist | UNCHECKED | 28.339 | 11.497 | 16.842 | 59.43% |

## Largest consistent reductions

| Scenario | CHECKED sum (ms) | UNCHECKED sum (ms) | Difference |
| --- | --- | --- | --- |
| interoperability.established_input.normalized_altitude_saliency | 76.517 | 55.865 | -26.99% |
| casf.established_input.incremental_sequence_int32_area | 1998.048 | 1574.739 | -21.19% |
| casf.established_input.incremental_area_heavy_step_after_medium | 50.582 | 40.729 | -19.48% |
| interoperability.established_input.topological_saliency | 105.532 | 86.212 | -18.31% |
| casf.established_input.incremental_sequence_float_bounding_box_diagonal | 1189.093 | 982.394 | -17.38% |
| casf.established_input.incremental_area_light_step | 134.929 | 111.930 | -17.05% |
| casf.established_input.incremental_area_sequence | 275.956 | 230.392 | -16.51% |
| casf.established_input.incremental_area_medium_step_after_light | 90.847 | 76.662 | -15.61% |
| filters.established_input.adaptive_mser_criterion | 5.547 | 4.756 | -14.25% |
| casf.established_input.incremental_bounding_box_diagonal_sequence | 304.538 | 261.874 | -14.01% |
| casf.end_to_end.pipeline_area_sequence | 376.628 | 323.990 | -13.98% |
| filters.established_input.extinction_ranked_formal_saliency | 437.539 | 381.031 | -12.92% |

## Largest increases

| Scenario | CHECKED sum (ms) | UNCHECKED sum (ms) | Difference |
| --- | --- | --- | --- |
| attribute_groups.established_input.moments | 16.246 | 17.862 | +9.95% |
| attributes.established_input.scalar_bitquads_length_average | 50.775 | 53.824 | +6.01% |
| filters.established_input.direct_criterion | 2.110 | 2.233 | +5.85% |
| attributes.established_input.scalar_bitquads_width_average | 51.125 | 53.864 | +5.36% |
| attributes.established_input.scalar_bitquads_perimeter_average | 50.867 | 53.564 | +5.30% |
| attributes.established_input.scalar_bitquads_perimeter | 50.679 | 53.328 | +5.23% |
| attributes.established_input.scalar_bitquads_number_holes | 50.692 | 53.274 | +5.09% |
| attribute_groups.established_input.boundary_sequential_scalars | 562.206 | 589.897 | +4.93% |

Small increases and reductions in sub-millisecond scalar scenarios should be treated as timing noise unless supported by a targeted experiment.

## Attribute grouping

| Mode | Group/bundle | Median scalars/grouped |
| --- | --- | --- |
| CHECKED | attribute_bundles.casf_area | 0.99x |
| CHECKED | attribute_bundles.casf_bounding_box | 2.34x |
| CHECKED | attribute_bundles.connected_filter | 1.30x |
| CHECKED | attribute_bundles.geometric_moments | 2.27x |
| CHECKED | attribute_bundles.max_distance_filter | 1.00x |
| CHECKED | attribute_bundles.moment_shape_filter | 5.87x |
| CHECKED | attribute_bundles.radiometric_size | 1.26x |
| CHECKED | attribute_bundles.shape_boundary | 2.59x |
| CHECKED | attribute_groups.all | 2.30x |
| CHECKED | attribute_groups.boundary | 7.86x |
| CHECKED | attribute_groups.gray_level | 2.46x |
| CHECKED | attribute_groups.moments | 6.68x |
| CHECKED | attribute_groups.shape | 2.33x |
| CHECKED | attribute_groups.tree_topology | 2.40x |
| UNCHECKED | attribute_bundles.casf_area | 1.00x |
| UNCHECKED | attribute_bundles.casf_bounding_box | 2.35x |
| UNCHECKED | attribute_bundles.connected_filter | 1.12x |
| UNCHECKED | attribute_bundles.geometric_moments | 2.27x |
| UNCHECKED | attribute_bundles.max_distance_filter | 1.00x |
| UNCHECKED | attribute_bundles.moment_shape_filter | 5.80x |
| UNCHECKED | attribute_bundles.radiometric_size | 1.27x |
| UNCHECKED | attribute_bundles.shape_boundary | 2.63x |
| UNCHECKED | attribute_groups.all | 2.28x |
| UNCHECKED | attribute_groups.boundary | 7.92x |
| UNCHECKED | attribute_groups.gray_level | 2.51x |
| UNCHECKED | attribute_groups.moments | 6.61x |
| UNCHECKED | attribute_groups.shape | 2.40x |
| UNCHECKED | attribute_groups.tree_topology | 2.44x |

## CASF

CASF aggregate timings and invariant edit counts for every workload are retained in `campaign-casf.csv`. The common-scenario table allows direct aggregation without mixing workload-specific bundles.

The complete campaign tables are `campaign-workloads.csv`, `campaign-scenarios.csv`, `campaign-common-scenarios.csv`, `campaign-scope-summary.csv`, `campaign-scope-pairs.csv`, `campaign-attribute-pairs.csv`, `campaign-attribute-summary.csv`, and `campaign-casf.csv`. Campaigns with captured repetitions also include `campaign-sample-statistics.csv`.

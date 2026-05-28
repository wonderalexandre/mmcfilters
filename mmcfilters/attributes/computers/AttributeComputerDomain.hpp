#pragma once

namespace mmcfilters::attributes::computers {

/**
 * @brief Stable execution domain declared by an attribute computer family.
 *
 * @details
 * The domain is part of the orchestration contract, not a mathematical
 * classification of the descriptor itself. `Topology` families can run from a
 * `MorphologicalTree` alone. `Altitude` families need a dense altitude span in
 * addition to topology.
 */
enum class AttributeComputerDomain {
    /// Computes from support, image geometry, or parent/child topology only.
    Topology,

    /// Computes from topology plus a typed altitude span.
    Altitude
};

} // namespace mmcfilters::attributes::computers

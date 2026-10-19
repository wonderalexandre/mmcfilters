#pragma once

#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <optional>

namespace mmcfilters {

template <AltitudeValue T> class ValuedMorphologicalTree;
template <class T> struct NativeHierarchyOwned;

namespace detail {

/**
 * @brief Selects the initial mutation token for a materialized hierarchy.
 *
 * Component trees and Higra imports use `numNodes + 1` to account for storage
 * initialization, non-root links, and finalization. The token is assigned once
 * after bulk materialization.
 */
enum class MaterializedVersionPolicy { CanonicalNative, PreserveLinkedConstruction };

template <AltitudeValue T>
[[nodiscard]] ValuedMorphologicalTree<T>
materializeNativeHierarchy(NativeHierarchyOwned<T>&& hierarchy,
                           std::optional<NodeId> preservedExternalNodeIdOffset = std::nullopt,
                           MaterializedVersionPolicy versionPolicy = MaterializedVersionPolicy::CanonicalNative);

/**
 * @brief Capability token for internal tree construction entry points.
 *
 * Several constructors must stay syntactically public because the library is
 * header-only and construction is assembled across multiple headers. They are
 * nevertheless not part of the user-facing construction API. Requiring this
 * tag makes that distinction explicit: the constructors can be called only by
 * code that is allowed to instantiate `MorphologicalTreeConstructionTag`.
 *
 * The tag carries no state and has no runtime cost. Its only purpose is access
 * control and documentation of construction ownership.
 */
class MorphologicalTreeConstructionTag {
  public:
    /**
     * @brief Copying is allowed so authorized construction code can forward the token.
     */
    constexpr MorphologicalTreeConstructionTag(const MorphologicalTreeConstructionTag&) noexcept = default;
    /**
     * @brief Assigns the supplied object state.
     *
     * @return Reference to the assigned object.
     */
    constexpr MorphologicalTreeConstructionTag& operator=(const MorphologicalTreeConstructionTag&) noexcept = default;

  private:
    /**
     * @brief Restricts tag creation to the internal native-buffer materializer.
     */
    constexpr MorphologicalTreeConstructionTag() noexcept = default;

    template <AltitudeValue T>
    friend ValuedMorphologicalTree<T> materializeNativeHierarchy(NativeHierarchyOwned<T>&&, std::optional<NodeId>, MaterializedVersionPolicy);
};

} // namespace detail

} // namespace mmcfilters

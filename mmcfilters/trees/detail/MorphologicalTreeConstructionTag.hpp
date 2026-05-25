#pragma once

namespace mmcfilters {

class MorphologicalTreeFactory;

namespace detail {

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
    constexpr MorphologicalTreeConstructionTag& operator=(const MorphologicalTreeConstructionTag&) noexcept = default;

private:
    /**
     * @brief Private default construction keeps tag creation inside the factory.
     */
    constexpr MorphologicalTreeConstructionTag() noexcept = default;

    friend class ::mmcfilters::MorphologicalTreeFactory;
};

} // namespace detail

} // namespace mmcfilters

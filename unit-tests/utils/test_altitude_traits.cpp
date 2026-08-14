#include "mmcfilters/utils/Altitude.hpp"
#include "mmcfilters/utils/Common.hpp"

#include <cassert>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

using namespace mmcfilters;

template <class T> void checkAltitudeSpanAlias() {
    static_assert(AltitudeValue<T>);
    static_assert(std::is_same_v<NodeAltitudeSpan<T>, std::span<const T>>);

    std::vector<T> values{T{0}, T{1}, T{2}};
    NodeAltitudeSpan<T> view(values);

    assert(view.size() == values.size());
    assert(view[0] == T{0});
    assert(view[1] == T{1});
    assert(view[2] == T{2});
}

int main() {
    static_assert(AltitudeValue<std::uint8_t>);
    static_assert(AltitudeValue<std::int32_t>);
    static_assert(AltitudeValue<float>);
    static_assert(AltitudeValue<double>);
    static_assert(!AltitudeValue<bool>);
    static_assert(!AltitudeValue<std::int64_t>);
    static_assert(!AltitudeValue<std::uint64_t>);
    static_assert(!AltitudeValue<std::vector<int>>);

    static_assert(std::is_same_v<AltitudeDifference<std::uint8_t>, std::int64_t>);
    static_assert(std::is_same_v<AltitudeDifference<std::int32_t>, std::int64_t>);
    static_assert(std::is_same_v<AltitudeDifference<float>, float>);
    static_assert(std::is_same_v<AltitudeDifference<double>, double>);
    static_assert(std::is_same_v<NodeAltitudeBuffer<std::int32_t>, std::vector<std::int32_t>>);
    static_assert(std::is_same_v<NodeAltitudeBuffer<std::uint8_t>, std::vector<std::uint8_t>>);

    checkAltitudeSpanAlias<std::uint8_t>();
    checkAltitudeSpanAlias<std::int32_t>();
    checkAltitudeSpanAlias<float>();
    checkAltitudeSpanAlias<double>();

    return 0;
}

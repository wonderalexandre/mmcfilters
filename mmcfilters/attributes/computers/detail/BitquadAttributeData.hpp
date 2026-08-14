#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail {

/// Bit position of the top-left sample in canonical row-major order.
inline constexpr std::size_t topLeftBit = 0;
/// Bit position of the top-right sample in canonical row-major order.
inline constexpr std::size_t topRightBit = 1;
/// Bit position of the bottom-left sample in canonical row-major order.
inline constexpr std::size_t bottomLeftBit = 2;
/// Bit position of the bottom-right sample in canonical row-major order.
inline constexpr std::size_t bottomRightBit = 3;

/** @brief Integer code of one canonical four-sample bitquad state. */
using BitquadCode = std::uint8_t;

/** @brief Spatially named bitquad state in canonical row-major order. */
struct BitquadState {
    bool topLeft = false;     ///< Foreground state of the top-left sample.
    bool topRight = false;    ///< Foreground state of the top-right sample.
    bool bottomLeft = false;  ///< Foreground state of the bottom-left sample.
    bool bottomRight = false; ///< Foreground state of the bottom-right sample.
};

/**
 * @brief Encodes a spatial state as the displayed word z3z2z1z0.
 * @param state Spatial bitquad state.
 * @return Code z0 + 2z1 + 4z2 + 8z3.
 */
[[nodiscard]] constexpr BitquadCode bitquadCode(BitquadState state) noexcept {
    return static_cast<BitquadCode>((state.topLeft ? BitquadCode{1} << topLeftBit : 0) | (state.topRight ? BitquadCode{1} << topRightBit : 0) |
                                    (state.bottomLeft ? BitquadCode{1} << bottomLeftBit : 0) | (state.bottomRight ? BitquadCode{1} << bottomRightBit : 0));
}

/** @brief Top-left coordinate of one cell in the framed bitquad grid. */
struct FramedBitquadCell {
    int topRow = 0;     ///< Row of the cell's top-left sample.
    int leftColumn = 0; ///< Column of the cell's top-left sample.

    /** @brief Compares two framed-cell coordinates. @return True when both coordinates are equal. */
    friend bool operator==(const FramedBitquadCell&, const FramedBitquadCell&) = default;
};

/** @brief One of the five nonempty bitquad families. */
enum class BitquadFamily : std::uint8_t { Q1, Q2, QD, Q3, Q4 };

/**
 * @brief Classifies a canonical bitquad code into a nonempty family.
 * @param code Canonical code in `[0, 15]`.
 * @return No value for code zero; otherwise its family.
 */
[[nodiscard]] constexpr std::optional<BitquadFamily> bitquadFamily(BitquadCode code) {
    switch (code) {
    case 0:
        return std::nullopt;
    case 1:
    case 2:
    case 4:
    case 8:
        return BitquadFamily::Q1;
    case 3:
    case 5:
    case 10:
    case 12:
        return BitquadFamily::Q2;
    case 6:
    case 9:
        return BitquadFamily::QD;
    case 7:
    case 11:
    case 13:
    case 14:
        return BitquadFamily::Q3;
    case 15:
        return BitquadFamily::Q4;
    default:
        throw std::out_of_range("Bitquad code must be in [0, 15].");
    }
}

/** @brief Signed five-family contribution introduced at one node before aggregation. */
struct BitquadFamilyIncrement {
    int q1 = 0; ///< Signed Q1 increment.
    int q2 = 0; ///< Signed Q2 increment.
    int qd = 0; ///< Signed QD increment.
    int q3 = 0; ///< Signed Q3 increment.
    int q4 = 0; ///< Signed Q4 increment.

    /** @brief Compares two family increments. @return True when every signed coordinate is equal. */
    friend bool operator==(const BitquadFamilyIncrement&, const BitquadFamilyIncrement&) = default;
};

/** @brief Aggregated counts of the five nonempty bitquad families. */
struct BitquadFamilyCounts {
    int q1 = 0; ///< Count of Q1 configurations.
    int q2 = 0; ///< Count of Q2 configurations.
    int qd = 0; ///< Count of QD configurations.
    int q3 = 0; ///< Count of Q3 configurations.
    int q4 = 0; ///< Count of Q4 configurations.

    /** @brief Compares two family-count values. @return True when every count is equal. */
    friend bool operator==(const BitquadFamilyCounts&, const BitquadFamilyCounts&) = default;
};

/** @brief Signed pre-aggregation histogram over canonical nonempty codes 1 through 15. */
struct NonemptyBitquadStateHistogramIncrement {
    std::array<int, 15> bins{}; ///< Bin `code - 1` stores the increment of canonical `code`.

    /** @brief Accesses one nonempty canonical code. @param code Code in `[1, 15]`. @return Mutable signed increment. */
    [[nodiscard]] int& count(BitquadCode code) {
        if (code == 0 || code > 15) {
            throw std::out_of_range("Nonempty bitquad increment code must be in [1, 15].");
        }
        return bins[static_cast<std::size_t>(code - 1)];
    }

    /** @brief Accesses one nonempty canonical code. @param code Code in `[1, 15]`. @return Signed increment. */
    [[nodiscard]] const int& count(BitquadCode code) const {
        if (code == 0 || code > 15) {
            throw std::out_of_range("Nonempty bitquad increment code must be in [1, 15].");
        }
        return bins[static_cast<std::size_t>(code - 1)];
    }

    /** @brief Compares two nonempty-state increments. @return True when every signed bin is equal. */
    friend bool operator==(const NonemptyBitquadStateHistogramIncrement&, const NonemptyBitquadStateHistogramIncrement&) = default;
};

/** @brief Aggregated histogram over canonical nonempty codes 1 through 15. */
struct NonemptyBitquadStateHistogram {
    std::array<int, 15> bins{}; ///< Bin `code - 1` stores the count of canonical `code`.

    /** @brief Accesses one nonempty canonical code. @param code Code in `[1, 15]`. @return Mutable count. */
    [[nodiscard]] int& count(BitquadCode code) {
        if (code == 0 || code > 15) {
            throw std::out_of_range("Nonempty bitquad histogram code must be in [1, 15].");
        }
        return bins[static_cast<std::size_t>(code - 1)];
    }

    /** @brief Accesses one nonempty canonical code. @param code Code in `[1, 15]`. @return Count. */
    [[nodiscard]] const int& count(BitquadCode code) const {
        if (code == 0 || code > 15) {
            throw std::out_of_range("Nonempty bitquad histogram code must be in [1, 15].");
        }
        return bins[static_cast<std::size_t>(code - 1)];
    }

    /** @brief Compares two aggregated nonempty histograms. @return True when every count is equal. */
    friend bool operator==(const NonemptyBitquadStateHistogram&, const NonemptyBitquadStateHistogram&) = default;
};

/** @brief Full aggregated histogram over canonical bitquad codes 0 through 15. */
struct BitquadStateHistogram {
    std::array<int, 16> bins{}; ///< Bin `code` stores the count of canonical `code`.

    /** @brief Accesses one canonical code. @param code Code in `[0, 15]`. @return Mutable count. */
    [[nodiscard]] int& count(BitquadCode code) {
        if (code > 15) {
            throw std::out_of_range("Bitquad histogram code must be in [0, 15].");
        }
        return bins[static_cast<std::size_t>(code)];
    }

    /** @brief Accesses one canonical code. @param code Code in `[0, 15]`. @return Count. */
    [[nodiscard]] const int& count(BitquadCode code) const {
        if (code > 15) {
            throw std::out_of_range("Bitquad histogram code must be in [0, 15].");
        }
        return bins[static_cast<std::size_t>(code)];
    }

    /** @brief Returns all bins in canonical code order. @return Read-only 16-bin span. */
    [[nodiscard]] std::span<const int, 16> values() const noexcept { return bins; }

    /** @brief Compares two full state histograms. @return True when every count is equal. */
    friend bool operator==(const BitquadStateHistogram&, const BitquadStateHistogram&) = default;
};

} // namespace mmcfilters::attributes::computers::detail

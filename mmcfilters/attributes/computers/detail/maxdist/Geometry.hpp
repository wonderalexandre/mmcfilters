#pragma once

#include <utility>
#include <ostream>

namespace mmcfilters::attributes::computers::detail::maxdist {
/**
 * @brief Integer point in image coordinates.
 *
 * The MAX_DIST implementation treats `x` as the column coordinate and `y`
 * as the row coordinate. This type is intentionally lightweight: it does
 * not know the image domain and therefore performs no bounds checking.
 * Bounds are enforced by Box2D when a point is converted to a linear pixel
 * index.
 */
class Point2D {
  public:
    /**
     * @brief Creates a point with explicit column (`x`) and row (`y`).
     *
     * @param x Horizontal coordinate.
     * @param y Vertical coordinate.
     */
    Point2D(int x, int y) : x_{x}, y_{y} {}

    /**
     * @brief Read-only coordinate accessors.
     *
     * @return Read-only coordinate accessors.
     */
    inline int x() const noexcept { return x_; }
    /**
     * @brief Returns or updates the vertical coordinate.
     *
     * @return Or updates the vertical coordinate.
     */
    inline int y() const noexcept { return y_; }

    /**
     * @brief Mutable coordinate accessors.
     *
     * These are used by low-level geometry code that updates temporary
     * points in place. The caller remains responsible for keeping the point
     * inside the intended image domain.
     *
     * @return Reference to the resulting object.
     */
    inline int& x() noexcept { return x_; }
    /**
     * @brief Returns or updates the vertical coordinate.
     *
     * @return Or updates the vertical coordinate.
     */
    inline int& y() noexcept { return y_; }

    /**
     * @brief Replaces one coordinate without touching the other.
     *
     * @param val Value used by the operation.
     */
    inline void x(int val) noexcept { x_ = val; }
    /**
     * @brief Returns or updates the vertical coordinate.
     *
     * @param val Value whose geometric representation is requested.
     */
    inline void y(int val) noexcept { y_ = val; }

    /**
     * @brief Returns component-wise point addition.
     *
     * @param other Object to compare with or transfer from.
     * @return Component-wise point addition.
     */
    Point2D add(const Point2D& other) const noexcept { return Point2D{x_ + other.x_, y_ + other.y_}; }

    /**
     * @brief Returns component-wise point subtraction.
     *
     * @param other Object to compare with or transfer from.
     * @return Component-wise point subtraction.
     */
    Point2D sub(const Point2D& other) const noexcept { return Point2D{x_ - other.x_, y_ - other.y_}; }

  private:
    /**
     * @brief Column coordinate.
     */
    int x_;

    /**
     * @brief Row coordinate.
     */
    int y_;
};

/**
 * @brief Component-wise point addition.
 *
 * @param p Point used by the operation.
 * @param q Second point used by the operation.
 * @return Component-wise sum of the two points.
 */
inline Point2D operator+(const Point2D& p, const Point2D& q) { return p.add(q); }

/**
 * @brief Component-wise point subtraction.
 *
 * @param p Point used by the operation.
 * @param q Second point used by the operation.
 * @return Component-wise difference of the two points.
 */
inline Point2D operator-(const Point2D& p, const Point2D& q) { return p.sub(q); }

/**
 * @brief Debug stream representation as `(x, y)`.
 *
 * @param os Output stream receiving the representation.
 * @param p Point used by the operation.
 * @return Reference to the resulting object.
 */
inline std::ostream& operator<<(std::ostream& os, const Point2D& p) { return os << "(" << p.x() << ", " << p.y() << ")"; }

/**
 * @brief Inclusive axis-aligned rectangular image domain.
 *
 * `topLeft` and `bottomRight` are both part of the domain. Linear indices
 * are row-major and always relative to `topLeft`, so shifted boxes satisfy
 * `index(point(i)) == i`. Constructors do not normalize invalid extents:
 * callers must pass `bottomRight >= topLeft` or positive width/height.
 */
class Box2D {
  public:
    /**
     * @brief Creates a box from inclusive corner points.
     *
     * @param topLeft Inclusive top-left domain corner.
     * @param bottomRight Inclusive bottom-right domain corner.
     */
    Box2D(const Point2D& topLeft, const Point2D& bottomRight)
        : topLeft_{topLeft}, bottomRight_{bottomRight}, width_{bottomRight.x() - topLeft.x() + 1}, height_{bottomRight.y() - topLeft.y() + 1} {}

    /**
     * @brief Creates a box from inclusive corner coordinates.
     *
     * @param topLeftX Horizontal coordinate of the top-left corner.
     * @param topLeftY Vertical coordinate of the top-left corner.
     * @param bottomRightX Horizontal coordinate of the bottom-right corner.
     * @param bottomRightY Vertical coordinate of the bottom-right corner.
     */
    Box2D(int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
        : topLeft_{topLeftX, topLeftY}, bottomRight_{bottomRightX, bottomRightY}, width_{bottomRightX - topLeftX + 1}, height_{bottomRightY - topLeftY + 1} {}

    /**
     * @brief Creates a box from its top-left corner and extent.
     *
     * @param topLeft Inclusive top-left domain corner.
     * @param width Width of the domain.
     * @param height Height of the domain.
     */
    Box2D(const Point2D& topLeft, int width, int height)
        : topLeft_{topLeft}, bottomRight_{topLeft_.x() + width - 1, topLeft_.y() + height - 1}, width_{width}, height_{height} {}

    /**
     * @brief Creates an origin-based image domain of size `width x height`.
     *
     * @param width Width of the domain.
     * @param height Height of the domain.
     */
    Box2D(int width, int height) : topLeft_{0, 0}, bottomRight_{width - 1, height - 1}, width_{width}, height_{height} {}

    /**
     * @brief Domain geometry accessors.
     *
     * @return Domain geometry accessors.
     */
    inline Point2D topLeft() const noexcept { return topLeft_; }
    /**
     * @brief Returns the inclusive bottom-right corner.
     *
     * @return The inclusive bottom-right corner.
     */
    inline Point2D bottomRight() const noexcept { return bottomRight_; }

    /**
     * @brief Returns the top-left horizontal coordinate.
     *
     * @return The top-left horizontal coordinate.
     */
    inline int topLeftX() const noexcept { return topLeft_.x(); }
    /**
     * @brief Returns the top-left vertical coordinate.
     *
     * @return The top-left vertical coordinate.
     */
    inline int topLeftY() const noexcept { return topLeft_.y(); }
    /**
     * @brief Returns the bottom-right horizontal coordinate.
     *
     * @return The bottom-right horizontal coordinate.
     */
    inline int bottomRightX() const noexcept { return bottomRight_.x(); }
    /**
     * @brief Returns the bottom-right vertical coordinate.
     *
     * @return The bottom-right vertical coordinate.
     */
    inline int bottomRightY() const noexcept { return bottomRight_.y(); }

    /**
     * @brief Returns the domain width.
     *
     * @return The domain width.
     */
    inline int width() const noexcept { return width_; }
    /**
     * @brief Returns the domain height.
     *
     * @return The domain height.
     */
    inline int height() const noexcept { return height_; }

    /**
     * @brief Converts a point to a row-major index relative to this box.
     *
     * @param p Point used by the operation.
     * @return The zero-based linear index, or `-1` when the point is outside
     * the inclusive box.
     *
     */
    int index(const Point2D& p) const noexcept {
        if (contains(p)) {
            Point2D q = p - topLeft_;
            return q.y() * width_ + q.x();
        }

        return -1;
    }

    /**
     * @brief Converts coordinates to a row-major index relative to this box.
     *
     * @param x Horizontal coordinate.
     * @param y Vertical coordinate.
     * @return The zero-based linear index, or `-1` when `(x, y)` is outside
     * the inclusive box.
     *
     */
    int index(int x, int y) const noexcept {
        if (contains(x, y)) {
            int qx = x - topLeft_.x();
            int qy = y - topLeft_.y();
            return qy * width_ + qx;
        }
        return -1;
    }

    /**
     * @brief Converts a row-major index into a point in this box.
     *
     * The method assumes `0 <= idx < width() * height()`. It is intentionally
     * unchecked because it is called inside tight image loops.
     *
     * @param idx Zero-based index used by the operation.
     * @return The converted row-major index into a point in this box.
     */
    Point2D point(int idx) const noexcept {
        int px = idx % width_ + topLeft_.x();
        int py = idx / width_ + topLeft_.y();

        return {px, py};
    }

    /**
     * @brief Same conversion as point(), returned as a standard pair.
     *
     * @param idx Zero-based index used by the operation.
     * @return Converted point coordinates.
     */
    std::pair<int, int> pointPair(int idx) const noexcept {
        int px = idx % width_ + topLeft_.x();
        int py = idx / width_ + topLeft_.y();

        return std::make_pair(px, py);
    }

    /**
     * @brief Checks whether a point lies inside the inclusive domain.
     *
     * @param p Point used by the operation.
     * @return True if a point lies inside the inclusive domain; otherwise false.
     */
    bool contains(const Point2D& p) const noexcept {
        if (topLeft_.x() <= p.x() && p.x() <= bottomRight_.x() && topLeft_.y() <= p.y() && p.y() <= bottomRight_.y())
            return true;

        return false;
    }

    /**
     * @brief Checks whether coordinates lie inside the inclusive domain.
     *
     * @param x Horizontal coordinate.
     * @param y Vertical coordinate.
     * @return True if coordinates lie inside the inclusive domain; otherwise false.
     */
    bool contains(int x, int y) const noexcept {
        if (topLeft_.x() <= x && x <= bottomRight_.x() && topLeft_.y() <= y && y <= bottomRight_.y())
            return true;

        return false;
    }

  private:
    /**
     * @brief Inclusive top-left corner.
     */
    Point2D topLeft_;

    /**
     * @brief Inclusive bottom-right corner.
     */
    Point2D bottomRight_;

    /**
     * @brief Number of columns in the box.
     */
    int width_;

    /**
     * @brief Number of rows in the box.
     */
    int height_;
};
} // namespace mmcfilters::attributes::computers::detail::maxdist

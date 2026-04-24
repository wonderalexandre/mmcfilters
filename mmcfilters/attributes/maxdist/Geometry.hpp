#pragma once 

#include <utility>
#include <ostream>

namespace mmcfilters
{
  namespace maxdist
  {
    class Point2D
    {
    public:
      Point2D(int x, int y)
        : x_{x}, y_{y}
      {}

      inline int x() const noexcept { return x_; }
      inline int y() const noexcept { return y_; }

      inline int &x() noexcept { return x_; }
      inline int &y() noexcept { return y_; }

      inline void x(int val) noexcept { x_ = val; }
      inline void y(int val) noexcept { y_ = val; }

      Point2D add(const Point2D &other) const noexcept {
        return Point2D{ x_ + other.x_, y_ + other.y_ };
      }

      Point2D sub(const Point2D &other) const noexcept {
        return Point2D { x_ - other.x_, y_ - other.y_ };
      }

    private:
      int x_;
      int y_;
    };

    inline Point2D operator+(const Point2D &p, const Point2D &q)
    {
      return p.add(q);
    }

    inline Point2D operator-(const Point2D &p, const Point2D &q)
    {
      return p.sub(q);
    }

    inline std::ostream &operator<<(std::ostream &os, const Point2D &p)
    {
      return os << "(" << p.x() << ", " << p.y() <<  ")";
    }


    class Box2D
    {
    public:
      Box2D(const Point2D &topLeft, const Point2D &bottomRight)
        : topLeft_{topLeft}, bottomRight_{bottomRight},
          width_{bottomRight.x() - topLeft.x() + 1},
          height_{bottomRight.y() - topLeft.y() + 1}
      {}

      Box2D(int topLeftX, int topLeftY, int bottomRightX, int bottomRightY)
        : topLeft_{topLeftX, topLeftY}, bottomRight_{bottomRightX, bottomRightY},
          width_{bottomRightX - topLeftX + 1},
          height_{bottomRightY - topLeftY + 1}
      {}

      Box2D(const Point2D &topLeft, int width, int height)
        :topLeft_{topLeft},
         bottomRight_{topLeft_.x() + width -1, topLeft_.y() + height - 1},
         width_{width},
         height_{height}
      {}

      Box2D(int width, int height)
        : topLeft_{0, 0}, bottomRight_{width - 1, height - 1},
          width_{width}, height_{height}
      {}

      inline Point2D topLeft() const noexcept { return topLeft_; }
      inline Point2D bottomRight() const noexcept { return bottomRight_; }

      inline int topLeftX() const noexcept { return topLeft_.x(); }
      inline int topLeftY() const noexcept { return topLeft_.y(); }
      inline int bottomRightX() const noexcept { return bottomRight_.x(); }
      inline int bottomRightY() const noexcept { return bottomRight_.y(); }
    
      inline int width() const noexcept { return width_; }
      inline int height() const noexcept { return height_; }

      int index(const Point2D &p) const noexcept {
        if (contains(p)) {
          Point2D q = p + topLeft_;
          return q.y() * width_ + q.x();
        }

        return -1;
      }

      int index(int x, int y) const noexcept {
        if (contains(x, y)) {
          int qx = topLeft_.x() + x;
          int qy = topLeft_.y() + y;
          return qy * width_ + qx;
        }
        return -1;
      }

      Point2D point(int idx) const noexcept {
        int px = idx % width_ + topLeft_.x();
        int py = idx / width_ + topLeft_.y();

        return {px, py};
      }

      std::pair<int, int> pointPair(int idx) const noexcept {
        int px = idx % width_ + topLeft_.x();
        int py = idx / width_ + topLeft_.y();

        return std::make_pair(px, py);
      }

      bool contains(const Point2D &p) const noexcept {
        if (topLeft_.x() <= p.x() && p.x() <= bottomRight_.x() &&
            topLeft_.y() <= p.y() && p.y() <= bottomRight_.y())
          return true;

        return false;
      }

      bool contains(int x, int y) const noexcept {
        if (topLeft_.x() <= x && x <= bottomRight_.x() &&
            topLeft_.y() <= y && y <= bottomRight_.y())
          return true;

        return false;
      }

    private:
      Point2D topLeft_;
      Point2D bottomRight_;
      int width_;
      int height_;
    };
  }
}

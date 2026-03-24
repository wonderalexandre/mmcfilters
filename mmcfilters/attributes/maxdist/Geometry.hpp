#pragma once 

#include <utility>

namespace mmcfilters
{
  namespace maxdist
  {
    class Point2D
    {
    public:
      Point2D(int x, int y);

      inline int x() const noexcept { return x_; }
      inline int y() const noexcept { return y_; }

      inline int &x() noexcept { return x_; }
      inline int &y() noexcept { return y_; }

      inline void x(int val) noexcept { x_ = val; }
      inline void y(int val) noexcept { y_ = val; }

      Point2D add(const Point2D &other) const noexcept;  
      Point2D sub(const Point2D &other) const noexcept;

    private:
      int x_;
      int y_;
    };

    Point2D operator+(const Point2D &p, const Point2D &q);
    Point2D operator-(const Point2D &p, const Point2D &q);

    class Box2D
    {
    public:
      Box2D(const Point2D &topLeft, const Point2D &bottomRight);
      Box2D(int topLeftX, int topLeftY, int bomttomRightX, int bototmRightY);
      Box2D(const Point2D &topLeft, int width, int height);
      Box2D(int width, int height);

      inline Point2D topLeft() const noexcept { return topLeft_; }
      inline Point2D bottomRight() const noexcept { return bottomRight_; }

      inline int topLeftX() const noexcept { return topLeft_.x(); }
      inline int topLeftY() const noexcept { return topLeft_.y(); }
      inline int bottomRightX() const noexcept { return bottomRight_.x(); }
      inline int bottomRightY() const noexcept { return bottomRight_.y(); }
    
      inline int width() const noexcept { return width_; }
      inline int height() const noexcept { return height_; }

      int index(const Point2D &p) const noexcept;
      int index(int x, int y) const noexcept;

      Point2D point(int idx) const noexcept;
      std::pair<int, int> pointPair(int idx) const noexcept;

      bool contains(const Point2D &p) const noexcept;
      bool contains(int x, int y) const noexcept;

    private:
      Point2D topLeft_;
      Point2D bottomRight_;
      int width_;
      int height_;
    };
  }
}
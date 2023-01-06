#ifndef __POINT_H__
#define __POINT_H__
#include <cmath>
                        // Class created instead of struct because of std::map find() method more efficient than in a std::unordered_map
                        // And we need to arrange our Points so that find() can work properly
class Point {
public:
    Point(double x, double y);
    ~Point();
    bool operator<(const Point& other) const;       // Overloading operators
    bool operator>(const Point& other) const;
    bool operator<=(const Point& other) const;
    bool operator>=(const Point& other) const;
    double getX() const;
    double getY() const;

private:
    double m_x, m_y;
};

#endif

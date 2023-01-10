#ifndef __POINT_H__
#define __POINT_H__
#include <cmath>

/**
* @class Point
* @brief Representation of a point in 2D space
*
* A class that defines a point in 2D space, with x and y coordinates.
* It also overloads various comparison operators for the point.
*/

class Point {
public:
    Point(double x, double y);
    ~Point();

    bool operator<(const Point& other) const;
    bool operator>(const Point& other) const;
    bool operator<=(const Point& other) const;
    bool operator>=(const Point& other) const;
    double getX() const;
    double getY() const;

private:
    double m_x, m_y;
};

#endif

#include "Point.h"

Point::Point(double x, double y) : m_x(x), m_y(y) {}

Point::~Point() {}

bool Point::operator<(const Point& other) const {                                                           // Operator overloading
    return std::sqrt(m_x * m_x + m_y * m_y) < std::sqrt(other.m_x * other.m_x + other.m_y * other.m_y);     // Distance to origin
}

bool Point::operator>(const Point& other) const {
    return other < *this;
}

bool Point::operator<=(const Point& other) const {
    return !(*this > other);
}

bool Point::operator>=(const Point& other) const {
    return !(*this < other);
}

double Point::getX() const {          // Getters for m_x and m_y
    return m_x;
}
double Point::getY() const {
    return m_y;
}
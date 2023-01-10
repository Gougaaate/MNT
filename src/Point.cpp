#include "Point.hpp"

Point::Point(double x, double y) : m_x(x), m_y(y) {}

Point::~Point() = default;          // Default destructor

/**
    * @brief Overloads the < operator for Point class
    * @param other Point object to compare with
    * @return true if the x and y coordinates of this point are less than the x and y coordinates of the other point
    *
    * The other functions do the same job with the others operators.
    * This overloading is used because o the std::map which needs to be sorted in order to compute in O(1) complexity.
    */

bool Point::operator<(const Point& other) const {
    return std::sqrt(m_x * m_x + m_y * m_y) < std::sqrt(other.m_x * other.m_x + other.m_y * other.m_y);     // Distance to origin used to compare the points among themselves
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
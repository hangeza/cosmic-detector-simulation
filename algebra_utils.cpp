#include <cmath>
#include <vector>

#include "algebra_types.h"
#include "algebra_utils.h"
#include "geometry_types.h"
#include <cassert>

std::ostream& operator<<(std::ostream& os, const std::valarray<double>& p)
{
    os << "(";
    for (double element : p) {
        os << element << " ";
    }
    os << ")";
    return os;
}

auto norm(const Vector& vec) -> double
{
    constexpr double exponent { 2 };
    return std::sqrt(std::pow(vec, exponent).sum());
}

auto cross_product(const Vector& a, const Vector& b) -> Vector
{
    assert(a.size() == 3);
    assert(b.size() == 3);
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0] };
}

bool inEpsilon(double value, double eps)
{
    return (std::fabs(value) <= eps);
}

bool isFuzzySame(const std::valarray<double>& a, const std::valarray<double>& b, double eps)
{
    assert(a.size() == b.size());
    if (!inEpsilon(norm(a - b), eps))
        return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!inEpsilon(a[i] - b[i]))
            return false;
    }
    return true;
}

double getBoundingBoxVolume(const std::pair<Point, Point>& bb)
{
    double volume {};
    double lx { bb.second[0] - bb.first[0] };
    double ly { bb.second[1] - bb.first[1] };
    double lz { bb.second[2] - bb.first[2] };
    volume = lx * ly * lz;
    return volume;
}

Point rotate(const Point& p, const Vector& rot_axis, double rot_angle)
{
    if (inEpsilon(rot_angle))
        return p;
    const Vector k { rot_axis / norm(rot_axis) };
    const double c { std::cos(rot_angle) };
    const double s { std::sin(rot_angle) };
    return { p + s * cross_product(k, p) + (1. - c) * cross_product(k, cross_product(k, p)) };
}

Vector operator*(const matrix2d<double>& lhs, const Vector& rhs)
{
    assert(lhs.cols() == rhs.size());
    Vector result(lhs.cols());
    for (size_t j = 0; j < lhs.cols(); ++j) {
        // Take dot product of row a[i] and col b[j]
        result[j] = (lhs.row(j) * rhs).sum();
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const matrix2d<double>& m)
{
    for (size_t r { 0 }; r < m.rows(); ++r) {
        os << "|";
        for (size_t c { 0 }; c < m.cols(); ++c) {
            os << m(r, c) << " ";
        }
        os << "|\n";
    }
    return os;
}

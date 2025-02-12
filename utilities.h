#pragma once

#include "algebra_utils.h"
#include <cmath>
#include <vector>

#if (__cplusplus < 202002L)
// mimic the identity function object for c++ versions prior to C++20 standard
namespace std {
template <class T>
struct identity {
    constexpr T&& operator()(T&& t) const noexcept { return std::forward<T>(t); }
    constexpr T operator()(const T& t) const noexcept { return t; }
};
}
#endif

template <typename T>
struct DataItem {
    T value {};
    T error {};
};

template <typename T, typename U>
using MeasurementVector = std::vector<std::pair<DataItem<T>, DataItem<U>>>;

constexpr double toDeg(double x) { return x * 180 / pi(); }
constexpr double toRad(double x) { return x * pi() / 180; }

void export_file(const MeasurementVector<double, double>& data, const std::string& filename);

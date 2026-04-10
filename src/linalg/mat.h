#include <array>

#pragma once

std::array<std::array<double, 3>, 3> inv3(const std::array<std::array<double, 3>, 3>& mat);
std::array<double, 3> mxv3(const std::array<std::array<double, 3>, 3>& m, const std::array<double, 3>& v);
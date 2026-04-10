#include <array>

#pragma once

std::array<double, 3> getgravaccel(const std::array<double, 3>& x, int id, double et);
std::array<double, 3> gettotalgravaccel(const std::array<double, 3>& x, double et);
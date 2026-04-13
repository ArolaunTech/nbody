#include <array>

#pragma once

void getstate(int id, double j2000, double* out);
void getlatlonmoonssb(double lat, double lon, double j2000, double* out);

void logstate(const std::array<double, 6>& state);
void logposearthref(const std::array<double, 3>& logpos, double et);
void logstateearthref(const std::array<double, 6>& logstate, double et);
#include <cmath>

#include "vector.h"

double dot3(const double* a, const double* b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double distance3(const double* a, const double* b) {
	return std::sqrt(
		(a[0] - b[0]) * (a[0] - b[0]) + 
		(a[1] - b[1]) * (a[1] - b[1]) + 
		(a[2] - b[2]) * (a[2] - b[2])
	);
}
#include <iostream>

#include LIB

#include "ephem.h"

void getstate(int id, double j2000, double* out) {
	double lt;

	spkez_c(id, j2000, "ECLIPJ2000", "NONE", 0, out, &lt);
}

void getlatlonmoonssb(double lat, double lon, double j2000, double* out) {
	double pos[3];
	double rotate[3][3];
	double state[6];

	srfrec_c(301, lon, lat, pos);
	pxform_c("MOON_ME", "ECLIPJ2000", j2000, rotate);
	mxv_c(rotate, pos, pos);

	getstate(301, j2000, state);

	for (int i = 0; i < 3; i++) {
		out[i] = state[i] + pos[i];
	}
}

void logstate(const std::array<double, 6>& state) {
	std::cout << "Position: <" << state[0] << ", " << state[1] << ", " << state[2] << "> km\n";
	std::cout << "Velocity: <" << state[3] << ", " << state[4] << ", " << state[5] << "> km/s\n";
}

void logposearthref(const std::array<double, 3>& logpos, double et) {
	double state[6];

	getstate(399, et, state);

	std::cout << "Position: <" << logpos[0] - state[0] << ", " << logpos[1] - state[1] << ", " << logpos[2] - state[2] << "> km\n";
}

void logstateearthref(const std::array<double, 6>& logstate, double et) {
	double state[6];

	getstate(399, et, state);

	std::cout << "Position: <" << logstate[0] - state[0] << ", " << logstate[1] - state[1] << ", " << logstate[2] - state[2] << "> km\n";
	std::cout << "Velocity: <" << logstate[3] - state[3] << ", " << logstate[4] - state[4] << ", " << logstate[5] - state[5] << "> km/s\n";
}
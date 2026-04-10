#include <unordered_map>
#include <cmath>

#include LIB

#include "gravity.h"
#include "../ephem/ephem.h"

std::array<double, 3> getgravaccel(const std::array<double, 3>& x, int id, double et) {
	const std::unordered_map<int, double> GMs = { // Horizons values, km^3 s^-2
		{ 10, 1.3271244004127942E+11},
		{  1, 2.2031868551400003E+04},
        {  2, 3.2485859200000000E+05},
        {301, 4.9028001184575496E+03},
        {399, 3.9860043550702266E+05},
        {  4, 4.2828375815756102E+04},
        {  5, 1.2671276409999998E+08},
        {  6, 3.7940584841799997E+07},
        {  7, 5.7945563999999985E+06},
        {  8, 6.8365271005803989E+06},
        {  9, 9.7550000000000000E+02}
	};

	double state[6];

	getstate(id, et, state);

	std::array<double, 3> dx = {x[0] - state[0], x[1] - state[1], x[2] - state[2]};

	double dist = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	double grav = GMs.at(id) / dist / dist;

	double fdir[3] = {0, 0, 0};

	// Add perturbations if Earth (effects are insignicant from other bodies)
	if (id == 399) {
		// Earth gravity model (J2 only)
		const double RE = 6378.135;
		const double J2 = 0.00108262545;
		double C = GMs.at(id) * J2 * RE * RE;

		double rotate[3][3];
		pxform_c("ITRF93", "ECLIPJ2000", et, rotate);

		double pole[3] = {rotate[0][2], rotate[1][2], rotate[2][2]};

		double z = pole[0] * dx[0] + pole[1] * dx[1] + pole[2] * dx[2];
		
		double xdir[3] = {dx[0] - z * pole[0], dx[1] - z * pole[1], dx[2] - z * pole[2]};

		double x = std::sqrt(xdir[0] * xdir[0] + xdir[1] * xdir[1] + xdir[2] * xdir[2]);

		xdir[0] /= x;
		xdir[1] /= x;
		xdir[2] /= x;

		for (int i = 0; i < 3; i++) {
			fdir[i] = 
				3 * z / std::pow(dist, 5) * pole[i] + 
				(1.5 / std::pow(dist, 5) - 7.5 * z * z / std::pow(dist, 7)) * dx[i];
		}

		fdir[0] *= -C;
		fdir[1] *= -C;
		fdir[2] *= -C;
	}

	return std::array<double, 3>{-dx[0] * grav / dist + fdir[0], -dx[1] * grav / dist + fdir[1], -dx[2] * grav / dist + fdir[2]};
}

std::array<double, 3> gettotalgravaccel(const std::array<double, 3>& x, double et) {
	const int NBODIES = 11;

	const int ids[NBODIES] = {
		10,  // Sun
		1,   // Mercury
		2,   // Venus
		301, // Moon
		399, // Earth
		4,   // Mars barycenter
		5,   // Jupiter barycenter
		6,   // Saturn barycenter
		7,   // Uranus barycenter
		8,   // Neptune barycenter
		9    // Pluto barycenter
	};

	double ax = 0, ay = 0, az = 0;
	for (int i = 0; i < NBODIES; i++) {
		std::array<double, 3> grav = getgravaccel(x, ids[i], et);

		ax += grav[0];
		ay += grav[1];
		az += grav[2];
	}

	return std::array<double, 3>{ax, ay, az};
}
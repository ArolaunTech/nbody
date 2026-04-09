#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <iomanip>
#include <limits>
#include <unordered_map>

#include LIB

#include "ephem/ephem.h"
#include "physics/de.h"
#include "linalg/mat.h"

std::array<double, 6> getmoonejecstate(double lat, double lon, double et, double heading, double ejecvel) {
	double pos[3];
	double r[3];
	double z[3] = {0, 0, 1};
	double N[3];
	double E[3];
	double dir[3];
	double rotate[3][3];
	double state[6];

	getlatlonmoonssb(lat * rpd_c(), lon * rpd_c(), et, pos);
	srfrec_c(301, lon * rpd_c(), lat * rpd_c(), r);
	ucrss_c(z, r, E);
	ucrss_c(r, E, N);

	for (int i = 0; i < 3; i++) {
		dir[i] = N[i] * std::cos(heading * rpd_c()) + E[i] * std::sin(heading * rpd_c());
	}

	pxform_c("MOON_ME", "ECLIPJ2000", et, rotate);
	mxv_c(rotate, dir, dir);

	getstate(301, et, state);

	return std::array<double, 6>{
		pos[0], pos[1], pos[2],
		state[3] + ejecvel * dir[0], state[4] + ejecvel * dir[1], state[5] + ejecvel * dir[2]
	};
}

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

std::vector<double> fsimple(const std::vector<double> & y) { // Depot differential equation, only considers gravity
	std::array<double, 3> a = gettotalgravaccel(std::array<double, 3>{y[1], y[2], y[3]}, y[0]);

	return std::vector<double>{1, y[4], y[5], y[6], a[0], a[1], a[2]};
}

std::vector<double> fprimer(const std::vector<double> & y) { // Rocket pod differential equation, also considers thrust
	/*======== y structure ========*/
	// 0     ET (seconds since J2000)
	// 1-3   POS (km)
	// 4-6   VEL (km)
	// 7-9   p
	// 10-12 p'
	//
	// Total: 13 elements

	const double delta = 1e-3;

	std::array<double, 3> ag = gettotalgravaccel(std::array<double, 3>{y[1], y[2], y[3]}, y[0]);
	std::array<double, 3> ag2 = gettotalgravaccel(std::array<double, 3>{y[1] + delta * y[7], y[2] + delta * y[8], y[3] + delta * y[9]}, y[0]);

	return std::vector<double>{
		// ET
		1,
		// POS
		y[4], y[5], y[6],
		// VEL
		ag[0], ag[1], ag[2],
		// p
		y[10], y[11], y[12],
		// p'
		(ag2[0] - ag[0]) / delta,
		(ag2[1] - ag[1]) / delta,
		(ag2[2] - ag[2]) / delta,
	};
}

std::array<double, 6> propagateNbody(const std::array<double, 6>& x0, double t0, double t1, int N, Method method) {
	std::vector<double> y0 = {
		t0,
		x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]
	};

	std::vector<double> yf = integrate(
		fsimple,
		y0,
		(t1 - t0) / N,
		N,
		method
	);

	return std::array<double, 6>{yf[1], yf[2], yf[3], yf[4], yf[5], yf[6]};
}

std::array<double, 6> stateattime(const std::vector<std::vector<double> >& orbit, double et) {
	std::size_t Npts = orbit.size();

	if (Npts == 1) {
		return std::array<double, 6>{
			orbit[0][1], orbit[0][2], orbit[0][3], 
			orbit[0][4], orbit[0][5], orbit[0][6]
		};
	}

	if (et < orbit[0][0]) {
		std::cerr << "Invalid time!\n";

		return std::array<double, 6>{
			orbit[0][1], orbit[0][2], orbit[0][3], 
			orbit[0][4], orbit[0][5], orbit[0][6]
		};
	}
	if (et > orbit[Npts - 1][0]) {
		std::cerr << "Invalid time!\n";

		return std::array<double, 6>{
			orbit[Npts - 1][1], orbit[Npts - 1][2], orbit[Npts - 1][3], 
			orbit[Npts - 1][4], orbit[Npts - 1][5], orbit[Npts - 1][6]
		};
	}

	std::size_t lo = 0;
	std::size_t hi = Npts - 2;

	while (hi > lo) {
		std::size_t mid = (lo + hi) / 2;

		if (orbit[mid + 1][0] < et) {
			lo = mid + 1;
		} else if (orbit[mid][0] > et) {
			hi = mid - 1;
		} else {
			lo = mid;
			break;
		}
	}

	std::vector<double> istate = orbit[lo];
	step(fsimple, istate, et - istate[0], METHOD_DOPRI8);

	return std::array<double, 6>{
		istate[1], istate[2], istate[3], 
		istate[4], istate[5], istate[6]
	};
}

double dot (double* a, double* b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double distance(double* a, double* b) {
	return std::sqrt(
		(a[0] - b[0]) * (a[0] - b[0]) + 
		(a[1] - b[1]) * (a[1] - b[1]) + 
		(a[2] - b[2]) * (a[2] - b[2])
	);
}

std::array<std::array<double, 6>, 6> estimatestatetransmatrix(
	double t0, double t1, const std::array<double, 6>& state0, int N
) {
	const double epsilon[6] = {1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6};

	std::array<std::array<double, 6>, 6> out;

	std::array<double, 6> state1 = propagateNbody(state0, t0, t1, N, METHOD_DOPRI8);
	std::array<double, 6> state0perturbed = state0;

	for (int i = 0; i < 6; i++) {
		state0perturbed[i] += epsilon[i];

		std::array<double, 6> state1perturbed = propagateNbody(state0perturbed, t0, t1, N, METHOD_DOPRI8);

		for (int j = 0; j < 6; j++) {
			out[j][i] = (state1perturbed[j] - state1[j]) / epsilon[i];
		}

		state0perturbed[i] -= epsilon[i];
	}

	return out;
}

int main() {
	// Config
	double inittime = 1000;
	double ejectime = 1000;

	double baselat = 0;
	double baselon = -45;

	double baseheading = 90;
	double ejecvel = 2.5;

	// Load
	std::cout << "Loading files...\n";
	
	furnsh_c("data/metakernel.tm");

	std::cout << "Files loaded!\n\n";

	// Simulate depot
	double state[6];
	getstate(399, inittime, state);

	const double dt = 60;

	std::vector<double> y0 = {inittime, state[0] + 7500, state[1], state[2], state[3], state[4] + 8.999, state[5] - 0.009};
	std::vector<std::vector<double> > ys = integraterecord(
		fsimple,
		y0,
		dt,
		15000,
		2,
		METHOD_DOPRI8
	);

	// Rocket pod
	std::array<double, 6> podstate = getmoonejecstate(baselat, baselon, ejectime, baseheading, ejecvel);

	std::vector<double> impulsetimes = {5000, 472000};

	std::array<double, 6> podburn0 = propagateNbody(podstate, ejectime, impulsetimes[0], 200, METHOD_DOPRI8);
	std::array<double, 6> podburn0init = podburn0;

	std::array<double, 6> depotburn1 = stateattime(ys, impulsetimes[1]);

	double mag;

	std::array<double, 6> podburn1;

	for (int it = 0; it < 50; it++) {
		podburn1 = propagateNbody(podburn0, impulsetimes[0], impulsetimes[1], 250, METHOD_DOPRI8);

		std::array<std::array<double, 6>, 6> phi = estimatestatetransmatrix(impulsetimes[0], impulsetimes[1], podburn0, 250);
		std::array<std::array<double, 3>, 3> phi12;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				phi12[i][j] = phi[i][j + 3];
			}
		}

		std::array<std::array<double, 3>, 3> phi12inv = inv3(phi12);
		std::array<double, 3> deltas;

		for (int i = 0; i < 3; i++) {
			deltas[i] = 
				phi12inv[i][0] * (depotburn1[0] - podburn1[0]) +
				phi12inv[i][1] * (depotburn1[1] - podburn1[1]) +
				phi12inv[i][2] * (depotburn1[2] - podburn1[2]);
		}

		mag = std::sqrt(deltas[0] * deltas[0] + deltas[1] * deltas[1] + deltas[2] * deltas[2]);
		double mult = mag < 0.01 ? 1 : 0.01 / mag;

		for (int i = 0; i < 3; i++) {
			podburn0[i + 3] += deltas[i] * mult;
		}

		if (mag < 1e-10) break;
	}

	podburn1 = propagateNbody(podburn0, impulsetimes[0], impulsetimes[1], 250, METHOD_DOPRI8);

	double dv1 = distance(&podburn0[3], &podburn0init[3]);
	double dv2 = distance(&podburn1[3], &depotburn1[3]);

	std::cout << "Burn 1: " << dv1 << " km/s\n";
	std::cout << "Burn 2: " << dv2 << " km/s\n";
	std::cout << "Total: " << dv1 + dv2 << " km/s\n";

	std::array<double, 3> impulse0, impulse1, prod, prod2;
	for (int i = 0; i < 3; i++) {
		impulse0[i] = (podburn0[i + 3] - podburn0init[i + 3]) / dv1;
		impulse1[i] = (depotburn1[i + 3] - podburn1[i + 3]) / dv2;
	}

	std::array<std::array<double, 6>, 6> phi = estimatestatetransmatrix(impulsetimes[0], impulsetimes[1], podburn0, 250);
	std::array<std::array<double, 3>, 3> phi11, phi12, phi21, phi22;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			phi11[i][j] = phi[i][j];
			phi12[i][j] = phi[i][j + 3];
			phi21[i][j] = phi[i + 3][j];
			phi22[i][j] = phi[i + 3][j + 3];
		}
	}

	std::array<std::array<double, 3>, 3> phi12inv = inv3(phi12);

	for (int i = 0; i < 3; i++) {
		prod[i] = 
			phi11[i][0] * impulse0[0] + 
			phi11[i][1] * impulse0[1] + 
			phi11[i][2] * impulse0[2];

		prod[i] = impulse1[i] - prod[i];
	}

	for (int i = 0; i < 3; i++) {
		prod2[i] = 
			phi12inv[i][0] * prod[0] + 
			phi12inv[i][1] * prod[1] + 
			phi12inv[i][2] * prod[2];
	}

	std::array<double, 3> dpburn1;

	for (int i = 0; i < 3; i++) {
		dpburn1[i] = 
			phi21[i][0] * impulse0[0] + 
			phi21[i][1] * impulse0[1] + 
			phi21[i][2] * impulse0[2] + 
			phi22[i][0] * prod2[0] +
			phi22[i][1] * prod2[1] +
			phi22[i][2] * prod2[2];
	}

	double dp0 = dot(impulse0.data(), prod2.data());
	double dp1 = dot(impulse1.data(), dpburn1.data());

	const double maxpderivative = 1e-5;

	if (dp0 > maxpderivative) dp0 = maxpderivative;
	if (dp0 < -maxpderivative) dp0 = -maxpderivative;
	if (dp1 > maxpderivative) dp1 = maxpderivative;
	if (dp1 < -maxpderivative) dp1 = -maxpderivative;

	std::cout << dp0 << " " << dp1 << "\n";

	//

	std::vector<double> y0pod = {
		impulsetimes[0],
		podburn0[0], podburn0[1], podburn0[2],
		podburn0[3], podburn0[4], podburn0[5],
		impulse0[0], impulse0[1], impulse0[2],
		prod2[0], prod2[1], prod2[2]
	};

	std::vector<std::vector<double> > orbito = integraterecord(
		fprimer,
		y0pod,
		dt,
		15000,
		2,
		METHOD_DOPRI8
	);

	// Write file
	std::cout << "Writing to file...\n";
	const std::vector<std::string> labels = {
		"time", "x", "y", "z", "vx", "vy", "vz"
	};

	std::ofstream output;

	output.open("res/output.txt");
	
	output << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (std::size_t i = 0; i < labels.size(); i++) {
		output << labels[i] << " ";
	}
	output << "\n";

	for (std::size_t i = 0; i < orbito.size(); i++) {
		for (std::size_t j = 0; j < orbito[i].size(); j++) {
			output << orbito[i][j] << " ";
		}
		output << "\n";
	}

	output.close();

	std::ofstream outputdepot;

	outputdepot.open("res/output2.txt");

	outputdepot << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (std::size_t i = 0; i < labels.size(); i++) {
		outputdepot << labels[i] << " ";
	}
	outputdepot << "\n";

	for (std::size_t i = 0; i < ys.size(); i++) {
		for (std::size_t j = 0; j < ys[i].size(); j++) {
			outputdepot << ys[i][j] << " ";
		}
		outputdepot << "\n";
	}

	outputdepot.close();
}
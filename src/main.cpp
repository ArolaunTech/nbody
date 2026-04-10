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
#include "physics/gravity.h"
#include "linalg/mat.h"
#include "linalg/vector.h"

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

std::array<double, 3> target(const std::array<double, 6>& state0, const std::array<double, 6>& state1, double t0, double t1) {
	std::array<double, 6> state0burn = state0;

	for (int it = 0; it < 50; it++) {
		std::array<double, 6> res = propagateNbody(state0burn, t0, t1, 250, METHOD_DOPRI8);

		std::array<std::array<double, 6>, 6> phi = estimatestatetransmatrix(t0, t1, state0burn, 250);
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
				phi12inv[i][0] * (state1[0] - res[0]) +
				phi12inv[i][1] * (state1[1] - res[1]) +
				phi12inv[i][2] * (state1[2] - res[2]);
		}

		double mag = std::sqrt(deltas[0] * deltas[0] + deltas[1] * deltas[1] + deltas[2] * deltas[2]);
		double mult = mag < 0.01 ? 1 : 0.01 / mag;

		for (int i = 0; i < 3; i++) {
			state0burn[i + 3] += deltas[i] * mult;
		}

		if (mag < 1e-10) break;
	}

	return std::array<double, 3>{state0burn[3] - state0[3], state0burn[4] - state0[4], state0burn[5] - state0[5]};
}

std::array<std::array<double, 3>, 2> targetprimer(
	const std::array<double, 3>& p0, 
	const std::array<double, 3>& p1,
	const std::array<double, 6>& state0,
	double t0,
	double t1
) {
	// Computes primer derivatives at start and end points.

	std::array<std::array<double, 6>, 6> phi = estimatestatetransmatrix(t0, t1, state0, 250);
	std::array<std::array<double, 3>, 3> phi11, phi12, phi21, phi22;
	std::array<double, 3> prod, prod2, dpburn1;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			phi11[i][j] = phi[i][j];
			phi12[i][j] = phi[i][j + 3];
			phi21[i][j] = phi[i + 3][j];
			phi22[i][j] = phi[i + 3][j + 3];
		}
	}

	std::array<std::array<double, 3>, 3> phi12inv = inv3(phi12);

	prod = mxv3(phi11, p0);

	for (int i = 0; i < 3; i++) {
		prod[i] = p1[i] - prod[i];
	}

	prod2 = mxv3(phi12inv, prod);

	for (int i = 0; i < 3; i++) {
		dpburn1[i] = 
			phi21[i][0] * p0[0] + 
			phi21[i][1] * p0[1] + 
			phi21[i][2] * p0[2] + 
			phi22[i][0] * prod2[0] +
			phi22[i][1] * prod2[1] +
			phi22[i][2] * prod2[2];
	}

	std::array<std::array<double, 3>, 2> out;

	out[0] = prod2;
	out[1] = dpburn1;

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

	// BFGS
	std::array<std::array<double, 2>, 2> invhessian;

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			invhessian[i][j] = (i == j) ? 1 : 0;
		}
	}

	std::array<double, 6> podstate = getmoonejecstate(baselat, baselon, ejectime, baseheading, ejecvel);
	std::vector<double> impulsetimes = {5000, 472000};
	std::array<double, 6> podburn0, podburn1;
	std::array<double, 3> impulse0, impulse1, pder0;

	double dp0 = 0, dp1 = 0;

	for (int it2 = 0; it2 < 10; it2++) {
		const double coeff = 2e5;

		impulsetimes[0] += dp0 * coeff;
		impulsetimes[1] += dp1 * coeff;

		std::array<double, 6> podburn0init = propagateNbody(podstate, ejectime, impulsetimes[0], 200, METHOD_DOPRI8);
		podburn0 = podburn0init;

		std::array<double, 6> depotburn1 = stateattime(ys, impulsetimes[1]);

		std::array<double, 3> burn = target(podburn0init, depotburn1, impulsetimes[0], impulsetimes[1]);

		for (int i = 0; i < 3; i++) {
			podburn0[i + 3] += burn[i];
		}

		podburn1 = propagateNbody(podburn0, impulsetimes[0], impulsetimes[1], 250, METHOD_DOPRI8);

		double dv1 = distance3(&podburn0[3], &podburn0init[3]);
		double dv2 = distance3(&podburn1[3], &depotburn1[3]);

		std::cout << impulsetimes[0] << " " << impulsetimes[1] << "\n";
		std::cout << "Burn 1: " << dv1 << " km/s\n";
		std::cout << "Burn 2: " << dv2 << " km/s\n";
		std::cout << "Total: " << dv1 + dv2 << " km/s\n";

		for (int i = 0; i < 3; i++) {
			impulse0[i] = (podburn0[i + 3] - podburn0init[i + 3]) / dv1;
			impulse1[i] = (depotburn1[i + 3] - podburn1[i + 3]) / dv2;
		}

		std::array<std::array<double, 3>, 2> pders = targetprimer(impulse0, impulse1, podburn0, impulsetimes[0], impulsetimes[1]);
		pder0 = pders[0];

		dp0 = dot3(impulse0.data(), pders[0].data());
		dp1 = dot3(impulse1.data(), pders[1].data());

		std::cout << dp0 << " " << dp1 << "\n";

		const double maxpderivative = 6e-4;

		if (dp0 > maxpderivative) dp0 = maxpderivative;
		if (dp0 < -maxpderivative) dp0 = -maxpderivative;
		if (dp1 > maxpderivative) dp1 = maxpderivative;
		if (dp1 < -maxpderivative) dp1 = -maxpderivative;
	}
	//

	std::vector<double> y0pod = {
		impulsetimes[0],
		podburn0[0], podburn0[1], podburn0[2],
		podburn0[3], podburn0[4], podburn0[5],
		impulse0[0], impulse0[1], impulse0[2],
		pder0[0], pder0[1], pder0[2]
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
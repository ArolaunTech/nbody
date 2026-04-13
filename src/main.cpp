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

std::array<double, 3> target(
	const std::array<double, 6>& state0, 
	const std::array<double, 6>& state1, 
	double t0, 
	double t1
) {
	std::array<double, 6> state0burn = state0;
	double err;

	for (int it = 0; it < 250; it++) {
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

		err = std::sqrt(
			(state1[0] - res[0]) * (state1[0] - res[0]) +
			(state1[1] - res[1]) * (state1[1] - res[1]) +
			(state1[2] - res[2]) * (state1[2] - res[2])
		);

		if (err < 1e-6) break;
	}

	if (err > 1e-6) std::cout << err << "target\n";

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

std::vector<double> multiimpulsegradients(
	const std::vector<double>& impulsetimes,
	const std::vector<std::array<double, 3> >& impulsepos,
	const std::array<double, 6>& podo,
	const std::array<double, 6>& depotf
) {
	std::size_t Nimpulses = impulsetimes.size();

	// Targetting and burns
	std::vector<std::array<std::array<double, 6>, 2> > legstates(Nimpulses - 1);
	std::vector<std::array<double, 3> > burns(Nimpulses);
	std::vector<std::array<double, 3> > primers(Nimpulses);

	// First leg
	std::array<double, 3> burn0 = target(
		podo, 
		std::array<double, 6>{
			impulsepos[1][0], 
			impulsepos[1][1], 
			impulsepos[1][2], 
			0, 0, 0
		}, 
		impulsetimes[0], 
		impulsetimes[1]
	);

	burns[0] = burn0;

	std::array<double, 6> curr = podo;

	curr[3] += burn0[0];
	curr[4] += burn0[1];
	curr[5] += burn0[2];

	legstates[0][0] = curr;
	legstates[0][1] = propagateNbody(legstates[0][0], impulsetimes[0], impulsetimes[1], 250, METHOD_DOPRI8);

	// Middle legs
	for (int i = 1; i < Nimpulses - 2; i++) {
		std::array<double, 3> burn = target(
			legstates[i - 1][1],
			std::array<double, 6>{
				impulsepos[i + 1][0], impulsepos[i + 1][1], impulsepos[i + 1][2],
				0, 0, 0
			},
			impulsetimes[i],
			impulsetimes[i + 1]
		);

		burns[i] = burn;

		for (int j = 0; j < 3; j++) {
			legstates[i][0][j] = impulsepos[i][j];
			legstates[i][0][j + 3] = legstates[i - 1][1][j + 3] + burn[j];
		}

		legstates[i][1] = propagateNbody(legstates[i][0], impulsetimes[i], impulsetimes[i + 1], 250, METHOD_DOPRI8);
	}

	// Last leg
	std::array<double, 3> burnlast = target(
		legstates[Nimpulses - 3][1],
		depotf,
		impulsetimes[Nimpulses - 2],
		impulsetimes[Nimpulses - 1]
	);

	burns[Nimpulses - 2] = burnlast;

	for (int j = 0; j < 3; j++) {
		legstates[Nimpulses - 2][0][j] = impulsepos[Nimpulses - 2][j];
		legstates[Nimpulses - 2][0][j + 3] = legstates[Nimpulses - 3][1][j + 3] + burnlast[j];
	}

	legstates[Nimpulses - 2][1] = propagateNbody(legstates[Nimpulses - 2][0], impulsetimes[Nimpulses - 2], impulsetimes[Nimpulses - 1], 250, METHOD_DOPRI8);
	
	// Last burn
	for (int i = 0; i < 3; i++) {
		burns[Nimpulses - 1][i] = depotf[i + 3] - legstates[Nimpulses - 2][1][i + 3];
	}

	double dv = 0;
	for (int i = 0; i < Nimpulses; i++) {
		double mag = std::sqrt(burns[i][0] * burns[i][0] + burns[i][1] * burns[i][1] + burns[i][2] * burns[i][2]);
		dv += mag;

		if (mag == 0) {
			primers[i][0] = 1;
			primers[i][1] = 0;
			primers[i][2] = 0;
		} else {
			for (int j = 0; j < 3; j++) {
				primers[i][j] = burns[i][j] / mag;
			}
		}
	}

	// Primer derivatives
	std::vector<std::array<std::array<double, 3>, 2> > pders(Nimpulses - 1);

	for (int i = 0; i < Nimpulses - 1; i++) {
		pders[i] = targetprimer(
			primers[i],
			primers[i + 1],
			legstates[i][0],
			impulsetimes[i],
			impulsetimes[i + 1]
		);
	}

	// derivatives of cost
	std::vector<std::array<double, 3> > dJdr(Nimpulses - 1);
	std::vector<double> dJdt(Nimpulses);

	dJdt[0] = -(pders[0][0][0] * burns[0][0] + pders[0][0][1] * burns[0][1] + pders[0][0][2] * burns[0][2]);
	dJdt[Nimpulses - 1] = -(
		pders[Nimpulses - 2][1][0] * burns[Nimpulses - 1][0] + 
		pders[Nimpulses - 2][1][1] * burns[Nimpulses - 1][1] + 
		pders[Nimpulses - 2][1][2] * burns[Nimpulses - 1][2]
	);

	for (int i = 1; i < Nimpulses - 1; i++) {
		dJdt[i] = 
			pders[i - 1][1][0] * legstates[i - 1][1][3] +
			pders[i - 1][1][1] * legstates[i - 1][1][4] +
			pders[i - 1][1][2] * legstates[i - 1][1][5] -
			pders[i][0][0] * legstates[i][0][3] -
			pders[i][0][1] * legstates[i][0][4] -
			pders[i][0][2] * legstates[i][0][5];

		dJdr[i][0] = pders[i][0][0] - pders[i-1][1][0];
		dJdr[i][1] = pders[i][0][1] - pders[i-1][1][1];
		dJdr[i][2] = pders[i][0][2] - pders[i-1][1][2];
	}

	// Output
	std::vector<double> out(4 * Nimpulses - 5);

	for (std::size_t i = 0; i < Nimpulses; i++) {
		out[i] = dJdt[i];
	}

	for (std::size_t i = 0; i < Nimpulses - 2; i++) {
		out[Nimpulses + 3 * i + 0] = dJdr[i + 1][0];
		out[Nimpulses + 3 * i + 1] = dJdr[i + 1][1];
		out[Nimpulses + 3 * i + 2] = dJdr[i + 1][2];
	}

	out[4 * Nimpulses - 6] = dv;

	return out;
}

std::vector<std::array<std::array<double, 6>, 2> > getlegstates(
	const std::vector<double>& impulsetimes,
	const std::vector<std::array<double, 3> >& impulsepos,
	const std::array<double, 6>& podo,
	const std::array<double, 6>& depotf
) {
	std::size_t Nimpulses = impulsetimes.size();

	// Targetting and burns
	std::vector<std::array<std::array<double, 6>, 2> > legstates(Nimpulses - 1);

	// First leg
	std::array<double, 3> burn0 = target(
		podo, 
		std::array<double, 6>{
			impulsepos[1][0], 
			impulsepos[1][1], 
			impulsepos[1][2], 
			0, 0, 0
		}, 
		impulsetimes[0], 
		impulsetimes[1]
	);

	std::array<double, 6> curr = podo;

	curr[3] += burn0[0];
	curr[4] += burn0[1];
	curr[5] += burn0[2];

	legstates[0][0] = curr;
	legstates[0][1] = propagateNbody(legstates[0][0], impulsetimes[0], impulsetimes[1], 250, METHOD_DOPRI8);

	// Middle legs
	for (int i = 1; i < Nimpulses - 2; i++) {
		std::array<double, 3> burn = target(
			legstates[i - 1][1],
			std::array<double, 6>{
				impulsepos[i + 1][0], impulsepos[i + 1][1], impulsepos[i + 1][2],
				0, 0, 0
			},
			impulsetimes[i],
			impulsetimes[i + 1]
		);

		for (int j = 0; j < 3; j++) {
			legstates[i][0][j] = impulsepos[i][j];
			legstates[i][0][j + 3] = legstates[i - 1][1][j + 3] + burn[j];
		}

		legstates[i][1] = propagateNbody(legstates[i][0], impulsetimes[i], impulsetimes[i + 1], 250, METHOD_DOPRI8);
	}

	// Last leg
	std::array<double, 3> burnlast = target(
		legstates[Nimpulses - 3][1],
		depotf,
		impulsetimes[Nimpulses - 2],
		impulsetimes[Nimpulses - 1]
	);

	for (int j = 0; j < 3; j++) {
		legstates[Nimpulses - 2][0][j] = impulsepos[Nimpulses - 2][j];
		legstates[Nimpulses - 2][0][j + 3] = legstates[Nimpulses - 3][1][j + 3] + burnlast[j];
	}

	legstates[Nimpulses - 2][1] = propagateNbody(legstates[Nimpulses - 2][0], impulsetimes[Nimpulses - 2], impulsetimes[Nimpulses - 1], 250, METHOD_DOPRI8);

	return legstates;
}

double norm(const std::vector<double>& x) {
	double sqnorm = 0;

	for (std::size_t i = 0; i < x.size(); i++) {
		sqnorm += x[i] * x[i];
	}

	return std::sqrt(sqnorm);
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
	// Constants
	std::array<double, 6> podstate = getmoonejecstate(baselat, baselon, ejectime, baseheading, ejecvel);

	// Multi-impulse optimization
	std::vector<double> impulsetimes = {5000, 100000, 472000};

	const std::size_t Nimpulses = impulsetimes.size();
	const std::size_t vars = 4 * Nimpulses - 6;

	// Get states
	std::array<double, 6> depotf = stateattime(ys, impulsetimes[Nimpulses - 1]);
	std::array<double, 6> podo = propagateNbody(podstate, ejectime, impulsetimes[0], 100, METHOD_DOPRI8);

	std::vector<std::array<double, 3> > impulsepos(Nimpulses); // Positions of middle impulses

	std::array<double, 3> burn0vel = target(podo, depotf, impulsetimes[0], impulsetimes[Nimpulses - 1]);
	std::array<double, 6> curr = podo;

	for (int i = 0; i < 3; i++) {
		curr[i + 3] += burn0vel[i];
	}

	for (int i = 1; i < Nimpulses - 1; i++) {
		std::array<double, 6> state = propagateNbody(curr, impulsetimes[0], impulsetimes[i], 250, METHOD_DOPRI8);

		for (int j = 0; j < 3; j++) {
			impulsepos[i][j] = state[j];
		}
	}

	for (int i = 0; i < Nimpulses; i++) {
		logposearthref(impulsepos[i], impulsetimes[i]);
	}

	// BFGS variables
	std::vector<std::vector<double> > H(vars, std::vector<double>(vars));

	for (std::size_t i = 0; i < vars; i++) {
		for (std::size_t j = 0; j < vars; j++) {
			H[i][j] = (i == j) ? 1 : 0;
		}
	}

	std::vector<double> grad = multiimpulsegradients(
		impulsetimes, 
		impulsepos,
		podo, depotf
	);

	// Iterations
	for (int it = 0; it < 1000; it++) {
		double sqnorm = 0;

		for (std::size_t i = 0; i < vars; i++) {
			sqnorm += grad[i] * grad[i];
		}

		if (sqnorm < 1e-12) break;

		// Search direction
		std::vector<double> p(vars);

		for (std::size_t i = 0; i < vars; i++) {
			p[i] = 0;

			for (std::size_t j = 0; j < vars; j++) {
				p[i] -= H[i][j] * grad[j];
			}
		}

		double gradTp = 0;

		for (std::size_t i = 0; i < vars; i++) {
			gradTp += grad[i] * p[i];
		}

		// Line search
		double a = 200 / norm(p);

		const double c1 = 1e-4;
		const double c2 = 0.9;

		double fx = grad[vars];

		std::vector<double> gradnew;
		double gradTpnew;

		std::vector<double> impulsetimesnew(Nimpulses);
		std::vector<std::array<double, 3> > impulseposnew(Nimpulses);

		double it2 = 0;
		double prev = 1e9;

		do {
			if (it2 > 20) break;
			it2++;

			a /= 2;

			for (std::size_t i = 0; i < Nimpulses; i++) {
				impulsetimesnew[i] = impulsetimes[i] + a * p[i];
			}

			for (std::size_t i = 1; i < Nimpulses - 1; i++) {
				impulseposnew[i][0] = impulsepos[i][0] + a * p[Nimpulses + 3 * i - 3 + 0];
				impulseposnew[i][1] = impulsepos[i][1] + a * p[Nimpulses + 3 * i - 3 + 1];
				impulseposnew[i][2] = impulsepos[i][2] + a * p[Nimpulses + 3 * i - 3 + 2];
			}

			std::array<double, 6> depotfnew = stateattime(ys, impulsetimesnew[Nimpulses - 1]);
			std::array<double, 6> podonew = propagateNbody(podstate, ejectime, impulsetimesnew[0], 100, METHOD_DOPRI8);

			gradnew = multiimpulsegradients(
				impulsetimesnew, 
				impulseposnew,
				podonew, depotfnew
			);

			gradTpnew = 0;

			for (std::size_t i = 0; i < vars; i++) {
				gradTpnew += gradnew[i] * p[i];
			}

			if (gradnew[vars] > prev) break;

			prev = gradnew[vars];

			std::cout << gradnew[vars] << " " << it2 << " " << it << "\n";
		} while (
			(gradnew[vars] >= fx + c1 * a * gradTp) ||
			(gradTpnew <= c2 * gradTp)
		);

		std::vector<double> s(vars);
		std::vector<double> y(vars);

		double sTy = 0;

		for (std::size_t i = 0; i < vars; i++) {
			s[i] = a * p[i];
			y[i] = gradnew[i] - grad[i];

			sTy += s[i] * y[i];
		}

		double r = 1 / sTy;

		std::vector<std::vector<double> > li(vars, std::vector<double>(vars));
		std::vector<std::vector<double> > ri(vars, std::vector<double>(vars));
		std::vector<std::vector<double> > prod(vars, std::vector<double>(vars));
		std::vector<std::vector<double> > prod2(vars, std::vector<double>(vars));

		for (std::size_t i = 0; i < vars; i++) {
			for (std::size_t j = 0; j < vars; j++) {
				li[i][j] = (i == j) ? 1 : 0;
				ri[i][j] = (i == j) ? 1 : 0;

				li[i][j] -= s[i] * y[j] * r;
				ri[i][j] -= y[i] * s[j] * r;
			}
		}

		for (std::size_t i = 0; i < vars; i++) {
			for (std::size_t j = 0; j < vars; j++) {
				prod[i][j] = 0;

				for (std::size_t k = 0; k < vars; k++) {
					prod[i][j] += li[i][k] * H[k][j];
				}
			}
		}

		for (std::size_t i = 0; i < vars; i++) {
			for (std::size_t j = 0; j < vars; j++) {
				prod2[i][j] = 0;

				for (std::size_t k = 0; k < vars; k++) {
					prod2[i][j] += prod[i][k] * ri[k][j];
				}
			}
		}

		for (std::size_t i = 0; i < vars; i++) {
			for (std::size_t j = 0; j < vars; j++) {
				H[i][j] = prod2[i][j] + s[i] * s[j] * r;
			}
		}

		grad = gradnew;

		impulsetimes = impulsetimesnew;
		impulsepos = impulseposnew;

		std::cout << grad[vars] << "\n";

		for (std::size_t i = 0; i < Nimpulses; i++) {
			std::cout << impulsetimes[i] << " ";
		}

		std::cout << "\n";

		std::array<double, 6> depotfend = stateattime(ys, impulsetimes[Nimpulses - 1]);
		std::array<double, 6> podoend = propagateNbody(podstate, ejectime, impulsetimes[0], 100, METHOD_DOPRI8);
		std::vector<std::array<std::array<double, 6>, 2> > legstates = getlegstates(impulsetimes, impulsepos, podoend, depotfend);

		logstateearthref(podoend, impulsetimes[0]);

		for (std::size_t i = 0; i < Nimpulses - 1; i++) {
			logstateearthref(legstates[i][0], impulsetimes[i]);
			std::cout << "---------------------\n";
			logstateearthref(legstates[i][1], impulsetimes[i + 1]);
		}

		logstateearthref(depotfend, impulsetimes[Nimpulses - 1]);

		std::cout << "\n";
	}

	// Write file
	std::cout << "Writing to file...\n";
	const std::vector<std::string> labels = {
		"time", "x", "y", "z", "vx", "vy", "vz"
	};

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
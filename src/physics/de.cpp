#include <iostream>

#include "de.h"

std::vector<double> integrate(
	const std::function<std::vector<double>(const std::vector<double>&)>& f, 
	std::vector<double> y0, 
	double dt, 
	int steps, 
	Method method
) {
	std::vector<double> y = y0;

	for (int i = 0; i < steps; i++) {
		step(f, y, dt, method);
	}

	return y;
}

std::vector<std::vector<double> > integraterecord(
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	std::vector<double> y0,
	double dt,
	int steps,
	int interval,
	Method method
) {
	std::vector<double> y = y0;
	std::vector<std::vector<double> > out;

	for (int i = 0; i < steps; i++) {
		if (i % interval == 0) {
			out.push_back(y);
		}
		
		if (i % 1000 == 0) {
			std::cout << i << "\n";
		}

		step(f, y, dt, method);
	}

	return out;
}

void step(
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	std::vector<double>& y,
	double dt,
	Method method
) {
	if (method == METHOD_EULER) {
		std::vector<double> dy = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += dy[i] * dt;
		}
	} else if (method == METHOD_MIDPOINT) {
		std::vector<double> k1, k2;

		k1 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += k1[i] * dt / 2;
		}

		k2 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += k2[i] * dt - k1[i] * dt / 2;
		}
	} else if (method == METHOD_RK4) {
		std::vector<double> k1, k2, k3, k4;

		k1 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += k1[i] * dt / 2;
		}

		k2 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += k2[i] * dt / 2 - k1[i] * dt / 2;
		}

		k3 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += k3[i] * dt - k2[i] * dt / 2;
		}

		k4 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt / 6 - k3[i] * dt;
		}
	} else if (method == METHOD_DOPRI5) {
		std::vector<double> k1, k2, k3, k4, k5, k6, y0;

		y0 = y;
		k1 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 0.2;
		}

		k2 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 3/40 + k2[i] * dt * 9/40;
		}

		k3 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 44/45 - k2[i] * dt * 56/15 + k3[i] * dt * 32/9;
		}

		k4 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 19372/6561 - k2[i] * dt * 25360/2187 + k3[i] * dt * 64448/6561 - k4[i] * dt * 212/729;
		}

		k5 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 9017/3168 - k2[i] * dt * 355/33 + k3[i] * dt * 46732/5247 + k4[i] * dt * 49/176 - k5[i] * dt * 5103/18656;
		}

		k6 = f(y);

		for (std::size_t i = 0; i < y.size(); i++) {
			y[i] = y0[i] + k1[i] * dt * 35/384 + k3[i] * dt * 500/1113 + k4[i] * dt * 125/192 - k5[i] * dt * 2187/6784 + k6[i] * dt * 11/84;
		}
	} else if (method == METHOD_DOPRI8) {
		const int steps = 13;

		const std::array<std::array<double, steps>, steps> tableau = {{
			{{}},
			{{1./18}},
			{{1./48, 1./16}},
			{{1./32, 0, 3./32}},
			{{5./16, 0, -75./64, 75./64}},
			{{3./80, 0, 0, 3./16, 3./20}},
			{{29443841./614563906, 0, 0, 77736538./692538347, -28693883./1125000000, 23124283./1800000000}},
			{{16016141./946692911, 0, 0, 61564180./158732637, 22789713./633445777, 545815736./2771057229, -180193667./1043307555}},
			{{39632708./573591083, 0, 0, -433636366./683701615, -421739975./2616292301, 100302831./723423059, 790204164./839813087, 800635310./3783071287}},
			{{246121993./1340847787, 0, 0, -37695042795./15268766246, -309121744./1061227803, -12992083./490766935, 6005943493./2108947869, 393006217./1396673457, 123872331./1001029789}},
			{{-1028468189./846180014, 0, 0, 8478235783./508512852, 1311729495./1432422823, -10304129995./1701304382, -48777925059./3047939560, 15336726248./1032824649, -45442868181./3398467696, 3065993473./597172653}},
			{{185892177./718116043, 0, 0, -3185094517./667107341, -477755414./1098053517, -703635378./230739211, 5731566787./1027545527, 5232866602./850066563, -4093664535./808688257, 3962137247./1805957418, 65686358./487910083}},
			{{403863854./491063109, 0, 0, -5068492393./434740067, -411421997./543043805, 652783627./914296604, 11173962825./925320556, -13158990841./6184727034, 3936647629./1978049680, -160528059./685178525, 248638103./1413531060}}
		}};

		const std::array<double, steps> b = {
			14005451./335480064,
			0,
			0,
			0,
			0,
			-59238493./1068277825,
			181606767./758867731,
			561292985./797845732,
			-1041891430./1371343529,
			760417239./1151165299,
			118820643./751138087,
			-528747749./2220607170,
			0.25
		};

		std::vector<double> y0;
		std::array<std::vector<double>, steps> ks;
		for (int i = 0; i < steps; i++) ks[i] = y;

		y0 = y;

		for (int i = 0; i < steps; i++) {
			y = y0;

			for (int j = 0; j < steps; j++) {
				for (std::size_t k = 0; k < y.size(); k++) {
					y[k] += ks[j][k] * dt * tableau[i][j];
				}
			}

			ks[i] = f(y);
		}

		y = y0;

		for (int i = 0; i < steps; i++) {
			for (std::size_t k = 0; k < y.size(); k++) {
				y[k] += ks[i][k] * dt * b[i];
			}
		}
	}
}
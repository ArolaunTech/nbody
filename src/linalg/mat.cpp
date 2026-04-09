#include "mat.h"

std::array<std::array<double, 3>, 3> inv3(const std::array<std::array<double, 3>, 3>& mat) {
	double determinant = 
		mat[0][0] * mat[1][1] * mat[2][2] +
		mat[0][1] * mat[1][2] * mat[2][0] +
		mat[0][2] * mat[1][0] * mat[2][1] - 
		mat[0][2] * mat[1][1] * mat[2][0] -
		mat[0][0] * mat[1][2] * mat[2][1] -
		mat[0][1] * mat[1][0] * mat[2][2];

	std::array<std::array<double, 3>, 3> cofactors;

	for (int i = 0; i < 3; i++) {
		int li = i == 0 ? 1 : 0;
		int hi = i == 2 ? 1 : 2;

		for (int j = 0; j < 3; j++) {
			int lj = j == 0 ? 1 : 0;
			int hj = j == 2 ? 1 : 2;

			cofactors[i][j] = mat[li][lj] * mat[hi][hj] - mat[li][hj] * mat[hi][lj];
		}
	}

	cofactors[0][1] *= -1;
	cofactors[1][0] *= -1;
	cofactors[1][2] *= -1;
	cofactors[2][1] *= -1;

	std::array<std::array<double, 3>, 3> inverse;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			inverse[j][i] = cofactors[i][j] / determinant;
		}
	}

	return inverse;
}
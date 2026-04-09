#include <vector>
#include <functional>

#pragma once

enum Method {
	METHOD_EULER, // Tends to gain energy
	METHOD_MIDPOINT, // Tends to lose energy
	METHOD_RK4, // Almost exact
	METHOD_DOPRI5, // Almost exact
	METHOD_DOPRI8
};

std::vector<double> integrate(
	const std::function<std::vector<double>(const std::vector<double>&)>& f, 
	std::vector<double> y0, 
	double dt, 
	int steps, 
	Method method
);

std::vector<std::vector<double> > integraterecord(
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	std::vector<double> y0,
	double dt,
	int steps,
	int interval,
	Method method
);

void step(
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	std::vector<double>& y,
	double dt,
	Method method
);
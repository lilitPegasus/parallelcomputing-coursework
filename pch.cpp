#include "pch.h"
#include <cmath>

double ApproximationObj::getApproximation(const double omegaLeft, const double omegaRight, const double omegaCenter, const double tauStep, const double hStep)
{
	return  tauStep * (omegaLeft - 2 * omegaCenter + omegaRight) / pow(hStep, 2) + tauStep * a * omegaCenter + tauStep * b * pow(omegaCenter, 1 / 2) + omegaCenter;
}

double ApproximationObj::getAccurate(double aX, double aT)
{
	double x = pow(sqrt(-b / a) + c * std::exp((8.0 * a * aT) / 12.0 + sqrt(4.0 * a / 12.0) * aX), - 0.5);
	return x;
}

double ApproximationObj::bottomLimit(double aX)
{
	return getAccurate(aX, 0);
}

double ApproximationObj::leftLimit(double aT)
{
	return getAccurate(0, aT);
}

double ApproximationObj::rightLimit(double aT)
{
	return getAccurate(1, aT);
}

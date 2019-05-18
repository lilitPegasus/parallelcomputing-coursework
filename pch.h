#ifndef PCH_H
#define PCH_H


class ApproximationObj
{
public:
	ApproximationObj(double paramA, double paramB, double paramC) :
		a(paramA), b(paramB), c(paramC) {}

	double getApproximation(const double omegaLeft, const double omegaRight, const double omegaCenter, const double tauStep, const double hStep);
	double getAccurate(double x, double t);
	double bottomLimit(double x);
	double leftLimit(double t);
	double rightLimit(double t);

private:
	const double a;
	const double b;
	const double c;
};

#endif

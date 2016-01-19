#include "cqm.hpp"

#include <iostream>
#include <iterator>

struct Kernel
{
	std::complex<double> operator()(std::complex<double> const &s)
	{
		return om_d / ((s+a)*(s+a) + om_d*om_d) * std::exp(-s*tau);
	}

	static const double om_d;
	static const double a;
	static const double tau;
};

const double Kernel::om_d = 100.0 * 2. * M_PI;
const double Kernel::a = 1.;
const double Kernel::tau = .2;


int main(void)
{
	unsigned N = 10000;
	double dt = 1e-4;
	double delta = 1e-12;

	std::vector<double> excitation(N, 0.0);
	excitation[0] = 1.0;

	CQM<2> cqm(N, dt, delta);
	cqm.eval(excitation.begin(), excitation.end(), Kernel());

	std::ostream_iterator<double> out(std::cout, "\n");
	std::copy(cqm.get_time_response().begin(), cqm.get_time_response().end(), out);

	return 0;
}


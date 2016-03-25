#include "cqm/cqm.hpp"

#include <iostream>
#include <iterator>

#include <Eigen/Dense>

typedef Eigen::Matrix<double, 2, 1> ExcType;
typedef Eigen::Matrix<double, 3, 1> RespType;

typedef Eigen::Matrix<std::complex<double>, 2, 1> LaplaceExcType;
typedef Eigen::Matrix<std::complex<double>, 3, 1> LaplaceRespType;

struct Kernel
{
	typedef Eigen::Matrix<std::complex<double>, 3, 2> Result;
	Result operator()(std::complex<double> const &s) const
	{
		Result M = Result::Zero();
		M(0,0) = M(1,0) = -(M(2,0) = om_d / ((s+a)*(s+a) + om_d*om_d) * std::exp(-s*tau));
		return M;
	}

	static const double om_d;
	static const double a;
	static const double tau;
};


const double Kernel::om_d = 100.0 * 2. * M_PI;
const double Kernel::a = 1.;
const double Kernel::tau = 0;


int main(void)
{
	unsigned N = 100;
	double dt = 1e-3;
	double delta = 1e-12;

	std::vector<ExcType> excitation(N, ExcType::Zero());
	excitation[0](0) = 1.0;
	excitation[0](1) = 2.0;

	CQM<2, ExcType, LaplaceExcType, RespType, LaplaceRespType> cqm(N, dt, delta);
	cqm.eval(excitation.begin(), excitation.end(), Kernel());

	std::ostream_iterator<RespType> out(std::cout, "\n\n");
	std::copy(cqm.get_time_response().begin(), cqm.get_time_response().end(), out);

	return 0;
}


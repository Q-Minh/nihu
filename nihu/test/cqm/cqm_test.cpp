#include <boost/math/constants/constants.hpp>

#include "nihu/cqm/cqm.hpp"

#include <iostream>
#include <iterator>

#include <Eigen/Core>

using namespace boost::math::double_constants;

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
//		M(0,0) = M(1,0) = -(M(2,0) = om_d / ((s+a)*(s+a) + om_d*om_d) * std::exp(-s*tau));
		M(0,0) = M(1,0) = -(M(2,0) = std::exp(-s*tau));
		return M;
	}

	static const double om_d;
	static const double a;
	static const double tau;
};


const double Kernel::om_d = 100.0 * two_pi;
const double Kernel::a = 1.;
const double Kernel::tau = 1e-2;


int main(void)
{
	unsigned N = 1000;
	double dt = 5e-4;
	double delta = 1e-12;

	std::vector<ExcType> excitation(N, ExcType::Zero());
	int k = 0;
	double f0 = 1e2;
	for (auto it = excitation.begin(); it != excitation.end(); ++it, ++k)
		(*it)(0) = sin(two_pi * f0 * k * dt);

	CQM<2, ExcType, LaplaceExcType, RespType, LaplaceRespType> cqm(N, dt, delta);
	cqm.eval(excitation.begin(), excitation.end(), Kernel());

	std::ostream_iterator<RespType> out(std::cout, "\n\n");
	std::copy(cqm.get_time_response().begin(), cqm.get_time_response().end(), out);

	return 0;
}


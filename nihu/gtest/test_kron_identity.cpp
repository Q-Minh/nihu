#include <gtest/gtest.h>

#include "nihu/fmm/kron_identity.hpp"
#include <Eigen/Core>

TEST(kron_identity, basic)
{
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
	static const size_t Dim = 3;

	Mat m(2, 2);
	m << 0, 1, 2, 3;
	Vec v(Dim * 2, 1);
	v << 1, 2, 3, 4, 5, 6;
	Vec res = NiHu::fmm::create_kron_identity<Dim>(m) * (2 * v);
	std::cout << res << std::endl;
}

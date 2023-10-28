/// \file laplace_3d_fmm.cpp
/// \brief implementation of class fmm::laplace_3d_fmm

#include "laplace_3d_fmm.hpp"
#include <boost/math/special_functions/binomial.hpp>

namespace NiHu
{
namespace fmm
{

void laplace_3d_cluster::set_expansion_length(size_t L)
{
	m_expansion_length = L;
}

size_t laplace_3d_cluster::get_expansion_length() const
{
	return m_expansion_length;
}

size_t laplace_3d_cluster::data_size() const
{
	return (m_expansion_length + 1) * (m_expansion_length + 1);
}

laplace_3d_cluster::multipole_t
laplace_3d_cluster::zero_multipole() const
{
	return multipole_t(data_size(), 1);
}

laplace_3d_cluster::local_t
laplace_3d_cluster::zero_local() const
{
	return local_t(data_size(), 1);
}

Eigen::Index
laplace_3d_cluster::linear_index(Eigen::Index n, Eigen::Index m) const
{
	if (n < 0 || m > n || -m > n)
		throw std::logic_error("Invalid index in hat matrix");
	return n * n + n + m;
}

laplace_3d_fmm::m2l::result_t
laplace_3d_fmm::m2l::operator()(cluster_t to, cluster_t from) const
{
	using boost::math::binomial_coefficient;

	result_t res(to.data_size(), from.data_size());

	typedef cluster_t::bounding_box_t::location_t x_t;
	x_t x = from.get_bounding_box().get_center() - to.get_bounding_box().get_center();
	double rho = x.norm();
	double phi = std::atan2(x(1), x(0));
	double theta = std::acos(x(2) / rho);

	std::complex<double> const J(0.0, 1.0);

	for (size_t j = 0; j <= to.get_expansion_length(); ++j)
	{
		for (int k = -int(j); k <= int(j); ++k)
		{
			Eigen::Index row = to.linear_index(j, k);

			for (size_t n = 0; n <= from.get_expansion_length(); ++n)
			{
				for (int m = -int(n); m <= int(n); ++m)
				{
					Eigen::Index col = from.linear_index(n, m);

					res(row, col) =
						std::pow(J, std::abs(k - m) - std::abs(k) - std::abs(m)) *
						Y(j + n, m - k, theta, phi) *
						std::pow(-1, n) *
						std::pow(rho, -int(j + n + 1)) *
						std::sqrt(binomial_coefficient<double>(unsigned(j + k + n - m), unsigned(j + k))) *
						std::sqrt(binomial_coefficient<double>(unsigned(n + m + j - k), unsigned(n + m)));
				}
			}
		}
	}

	return res;
}

laplace_3d_fmm::m2m::result_t
laplace_3d_fmm::m2m::operator()(cluster_t to, cluster_t from) const
{
	using boost::math::binomial_coefficient;
	result_t res = result_t::Zero(to.data_size(), from.data_size());

	typedef cluster_t::bounding_box_t::location_t x_t;
	x_t x = from.get_bounding_box().get_center() - to.get_bounding_box().get_center();
	std::complex<double> const I(0.0, 1.0);

	double r = x.norm();
	double phi = std::atan2(x(1), x(0));
	double theta = std::acos(x(2) / r);

	int p = int(to.get_expansion_length());

	for (int j = 0; j <= p; ++j)
	{
		for (int k = -j; k <= +j; ++k)
		{
			Eigen::Index row_idx = to.linear_index(j, k);

			for (int n = 0; n <= j; ++n)
			{
				for (int m = -n; m <= +n; ++m)
				{
					if (n - m > j - k || m + n > j + k)
						continue;
					Eigen::Index col_idx = from.linear_index(j - n, k - m);

					res(row_idx, col_idx) =
						std::pow(I, std::abs(k) - std::abs(m) - std::abs(k - m)) *
						std::pow(r, n) *
						Y(n, -m, theta, phi) *
						std::sqrt(
							binomial_coefficient<double>(j - k, n - m) *
							binomial_coefficient<double>(j + k, n + m)
						);
				}
			}
		}
	}

	return res;
}

laplace_3d_fmm::l2l::result_t
laplace_3d_fmm::l2l::operator()(cluster_t to, cluster_t from) const
{
	using boost::math::binomial_coefficient;
	result_t res = result_t::Zero(to.data_size(), from.data_size());

	typedef cluster_t::bounding_box_t::location_t x_t;
	x_t x = from.get_bounding_box().get_center() - to.get_bounding_box().get_center();
	std::complex<double> const I(0.0, 1.0);
	double r = x.norm();
	double phi = std::atan2(x(1), x(0));
	double theta = std::acos(x(2) / r);

	int p = int(to.get_expansion_length());

	for (int j = 0; j <= p; ++j)
	{
		for (int k = -j; k <= +j; ++k)
		{
			Eigen::Index row_idx = to.linear_index(j, k);
			for (int n = j; n <= p; ++n)
			{
				for (int m = -n; m <= +n; ++m)
				{
					if (n - m < j - k || n + m < j + k)
						continue;
					Eigen::Index col_idx = from.linear_index(n, m);
					res(row_idx, col_idx) =
						std::pow(I, std::abs(m) - std::abs(m - k) - std::abs(k)) *
						Y(n - j, m - k, theta, phi) *
						std::pow(r, n - j) *
						std::pow(-1, n + j) *
						std::sqrt(
							binomial_coefficient<double>(n - m, j - k) *
							binomial_coefficient<double>(n + m, j + k)
						);
				}
			}
		}
	}

	return res;
}

laplace_3d_fmm::m2m
laplace_3d_fmm::create_m2m() const
{
	return m2m();
}

laplace_3d_fmm::l2l
laplace_3d_fmm::create_l2l() const
{
	return l2l();
}

laplace_3d_fmm::m2l
laplace_3d_fmm::create_m2l() const
{
	return m2l();
}


} // end of namespace fmm
} // namespace NiHu

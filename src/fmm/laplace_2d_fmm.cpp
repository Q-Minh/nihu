#include "laplace_2d_fmm.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace NiHu
{
namespace fmm
{



std::complex<double> laplace_2d_fmm::center2complex(cluster_t const& c)
{
	auto const& x = c.get_bounding_box().get_center();
	return std::complex<double>(x(0), x(1));
}

laplace_2d_fmm::m2m
laplace_2d_fmm::create_m2m() const
{
	return m2m();
}

laplace_2d_fmm::l2l
laplace_2d_fmm::create_l2l() const
{
	return l2l();
}

laplace_2d_fmm::m2l
laplace_2d_fmm::create_m2l() const
{
	return m2l();
}

void laplace_2d_cluster::set_expansion_length(size_t L)
{
	m_expansion_length = L;
}

size_t laplace_2d_cluster::get_expansion_length() const
{
	return m_expansion_length;
}

laplace_2d_cluster::multipole_t laplace_2d_cluster::zero_multipole() const
{
	return multipole_t::Zero(m_expansion_length + 1, 1);
}

laplace_2d_cluster::local_t laplace_2d_cluster::zero_local() const
{
	return local_t::Zero(m_expansion_length + 1, 1);
}

laplace_2d_fmm::m2l::result_t 
laplace_2d_fmm::m2l::operator()(cluster_t to, cluster_t from) const
{
	std::complex<double> z0 = center2complex(from) - center2complex(to);
	result_t res(to.get_expansion_length() + 1, from.get_expansion_length() + 1);
	res(0, 0) = std::log(-z0);
	for (size_t k = 1; k <= from.get_expansion_length(); ++k)
		res(0, k) = 1. / std::pow(-z0, k);
	for (size_t l = 1; l <= to.get_expansion_length(); ++l)
	{
		res(l, 0) = -1. / l / std::pow(z0, l);
		for (size_t k = 1; k <= from.get_expansion_length(); ++k)
			res(l, k) = 1. / std::pow(z0, l + k) * std::pow(-1, k)
			* boost::math::binomial_coefficient<double>(unsigned(k + l - 1), unsigned(k - 1));
	}
	return res;
}

laplace_2d_fmm::m2m::result_t 
laplace_2d_fmm::m2m::operator()(cluster_t to, cluster_t from) const
{
	std::complex<double> z0 = center2complex(from) - center2complex(to);

	result_t res = result_t::Zero(to.get_expansion_length() + 1, from.get_expansion_length() + 1);
	res(0, 0) = 1;
	for (size_t l = 1; l <= to.get_expansion_length(); ++l)
	{
		res(l, 0) = -std::pow(z0, l) / double(l);
		for (size_t k = 1; k <= l; ++k)
			res(l, k) = std::pow(z0, l - k) 
			* boost::math::binomial_coefficient<double>(unsigned(l - 1), unsigned(k - 1));
	}
	return res;
}

laplace_2d_fmm::l2l::result_t 
laplace_2d_fmm::l2l::operator()(cluster_t to, cluster_t from) const
{
	std::complex<double> z0 = center2complex(from) - center2complex(to);
	result_t res = result_t::Zero(to.get_expansion_length() + 1, from.get_expansion_length() + 1);
	for (size_t l = 0; l <= to.get_expansion_length(); ++l)
	{
		for (size_t k = l; k <= from.get_expansion_length(); ++k)
			res(l, k) = boost::math::binomial_coefficient<double>(unsigned(k), unsigned(l))
			* std::pow(-z0, k - l);
	}
	return res;
}

} // end of namespace fmm
} // namespace NiHu
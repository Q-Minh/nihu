#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <iostream>
#include "domain.hpp"

template <class Domain, unsigned N>
class GaussQuad;

template <unsigned N>
class GaussQuad<line_domain, N>
{
public:
	typedef Eigen::Matrix<double, N, 1> xi_type;
	typedef Eigen::Matrix<double, N, 1> weight_type;

	static void init(void)
	{
		typedef Eigen::Matrix<double, N, N> Mat_t;
		Mat_t A = Mat_t::Zero();
		for (unsigned i = 1; i < N; ++i)
			A(i, i-1) = A(i-1, i) = i / sqrt(4.0*(i*i)-1.0);
		Eigen::SelfAdjointEigenSolver<Mat_t> S(A);
		xi = S.eigenvalues();
		w = 2.0 * S.eigenvectors().row(0).cwiseAbs2().transpose();
	}
	
	static xi_type const &get_xi(void) { return xi; }
	static weight_type const &get_weight(void) { return w; }

protected:
	static xi_type xi;
	static weight_type w;
};

template <unsigned N>
typename GaussQuad<line_domain, N>::xi_type GaussQuad<line_domain, N>::xi;
template <unsigned N>
typename GaussQuad<line_domain, N>::weight_type GaussQuad<line_domain, N>::w;

int main(void)
{
	typedef GaussQuad<line_domain, 10> G;
	G::init();

	std::cout << G::get_xi() << std::endl << std::endl << G::get_weight() << std::endl;
	
	return 0;
}


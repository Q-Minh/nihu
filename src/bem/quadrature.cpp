#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <iostream>

#include "domain.hpp"

template <class Domain>
class GaussQuad;

template <>
class GaussQuad<line_domain>
{
public:
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec_t;

	GaussQuad(unsigned n)
	{
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat_t;
		Mat_t A = Mat_t::Zero(n,n);
		for (unsigned i = 1; i < n; ++i)
			A(i, i-1) = A(i-1, i) = i / sqrt(4.0*(i*i)-1.0);
		Eigen::SelfAdjointEigenSolver<Mat_t> S(A);
		xi = S.eigenvalues();
		w = 2.0 * S.eigenvectors().row(0).cwiseAbs2().transpose();
	}
	
	Vec_t const &get_xi(void) const { return xi; }
	Vec_t const &get_weight(void) const { return w; }

protected:
	Vec_t xi, w;
};

int main(void)
{
	GaussQuad<line_domain> G(10);
	std::cout << G.get_xi() << std::endl << std::endl << G.get_weight() << std::endl;
	
	return 0;
}


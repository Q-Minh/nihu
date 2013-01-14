#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
using namespace Eigen;

#include <iostream>

int main(void)
{
	typedef Matrix<double, Dynamic, Dynamic> Mat_t;
	typedef internal::traits<Mat_t >::Index Index;

	unsigned n = 5;
	Mat_t A = Mat_t::Zero(n,n);
	for (Index i = 1; i < n; ++i)
		A(i, i-1) = A(i-1, i) = i / sqrt(4.0*i*i - 1.0);
	SelfAdjointEigenSolver<Mat_t > S(A);
	Mat_t v = S.eigenvectors();
	std::cout << 2.0 * v.row(1).cwiseAbs2() << std::endl;
	
	return 0;
	
}


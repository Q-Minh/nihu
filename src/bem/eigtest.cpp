#include <Eigen/Dense>

typedef Eigen::Matrix<double, 100, 100> mat_t;

__forceinline mat_t add(mat_t const &a, mat_t const &b)
{
	return a + b;
}

#include <iostream>
using std::cout;


int main(void)
{
	mat_t a = mat_t::Random();
	mat_t b = mat_t::Random();
	auto c = a + b;
	cout << c(0,0);
	return 0;
}


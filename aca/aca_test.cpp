#include <iostream>
#include <complex>

#include "aca.hpp"

struct MAT
{
	typedef double Scalar;

	Scalar operator()(int i, int j) const 
	{
		return Scalar(i+j+100.);
	}
};

struct Laplace
{
	typedef double Scalar;

	Scalar operator()(int i, int j) const 
	{
		return 1./(i/10. - j/10.+2.);
	}
};

int main(void)
{
	int I[] = {0, 1, 2, 3, 4};
	int J[] = {5, 6, 7, 8};

	auto block = createBlock(Laplace(), &I[0], &J[0]);

	ACA<decltype(block), 5> aca(block, 5, 4);

	aca.process(1e-3);
	std::cout << aca.get_U() * aca.get_V().transpose() << '\n';

	std::cout << "rank: " << aca.get_U().cols() << '\n';

	std::cout << aca.get_U() << "\n\n" << aca.get_V() << std::endl;

	return 0;
}


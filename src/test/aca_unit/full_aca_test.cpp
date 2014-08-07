#include <iostream>
#include "aca/aca.hpp"

class Laplace
{
public:
	typedef double Scalar;

	Laplace(double x0, double Lx, int Nx, double y0, double Ly, int Ny)
		: x0(x0), dx(Lx/Nx), y0(y0), dy(Ly/Ny)
	{
	}

	double operator()(int i, int j) const
	{
		double x = x0 + i * dx;
		double y = y0 + j * dy;
		return 1./std::abs(y-x);
	}
private:
	double x0, dx, y0, dy;
};

int const Nx=200, Ny=100;
#define MULTILEVEL

int main(void)
{
	Laplace M(0., 1., Nx, 2., 1., Ny);

#ifdef MULTILEVEL
	Eigen::Matrix<int, 2, 2> rowClusters;
	rowClusters << 0,Nx/2,  Nx/2,Nx/2;

	Eigen::Matrix<int, 1, 2> colClusters;
	colClusters << 0,Ny;

	Eigen::Matrix<int, 2, 2> blocks;
	blocks << 0,0,  1,0;
#else
	Eigen::Matrix<int, 1, 2> rowClusters;
	rowClusters << 0,Nx;

	Eigen::Matrix<int, 1, 2> colClusters;
	colClusters << 0,Ny;

	Eigen::Matrix<int, 1, 2> blocks;
	blocks << 0,0;
#endif

	Eigen::VectorXd input(Ny,1), output(Nx,1);
	input.setRandom();
	output.setZero();

	ACA aca;
	aca.multiply(M, rowClusters, colClusters, blocks, input, output, 1e-3, 50);

	Eigen::VectorXd output_full(Nx,1);
	output_full.setZero();
	for (int i = 0; i < Nx; ++i)
		for (int j = 0; j < Ny; ++j)
			output_full(i) += M(i,j) * input(j);

	std::cout << "#elements: " << aca.get_nEval() << " / " <<
	aca.get_matrixSize() << " : " << aca.get_sizeCompression() << std::endl;
	std::cout << "error: " << (output - output_full).norm()/output.norm() << std::endl;

	return 0;
}


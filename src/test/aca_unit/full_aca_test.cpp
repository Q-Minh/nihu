#include <iostream>
#include "aca/aca.hpp"

int main(void)
{
	Eigen::MatrixXd M;
	M.resize(100, 100);
	for (int i = 0; i < 100; ++i)
		for (int j = 0; j < 100; ++j)
			M(i,j) = std::abs(1./(i+150-j));
			
	Eigen::Matrix<int, 6, 2> rowClusters(6,2);
	rowClusters << 50,50, 0,50, 25,25, 0,25, 75,25, 50,25;

	Eigen::Matrix<int, 6, 2> colClusters;
	colClusters << 0,50, 50,50, 0,25, 25,25, 50,25, 75,25;
	
	Eigen::Matrix<int, 6, 2> blocks(6,2);
	blocks << 0,0, 1,1, 2,2, 3,3, 4,4, 5,5;
	
	Eigen::VectorXd input(100,1), output(100,1);
	input.setRandom();
	output.setZero();
	
	ACA aca;
	aca.multiply(M, rowClusters, colClusters, blocks, input, output, 1e-5, 50);
		
	auto output_full = (M * input).eval();

	std::cout << "#elements: " << aca.get_nEval() << " / " <<
	aca.get_matrixSize() << " : " << aca.get_sizeCompression() << std::endl;
	std::cout << "error: " << (output - output_full).norm()/output.norm() << std::endl;

	return 0;
}


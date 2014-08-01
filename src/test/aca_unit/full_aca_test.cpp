#include <iostream>
#include <iterator>

#include "aca/aca.hpp"

int main(void)
{
	Eigen::MatrixXd M;
	M.resize(100, 100);
	for (int i = 0; i < 100; ++i)
		for (int j = 0; j < 100; ++j)
			M(i,j) = std::abs(1./(i+200-j));
			
	Eigen::MatrixXd U(100,50), V(100,50);
	
	double eps = 1e-7;

	int r = low_rank_approx(M, M.rows(), M.cols(), eps, 50, U, V);
	
	std::cout << "approximate rank: " << r << '\n';
	std::cout << "prescribed error: " << eps << '\n';
	std::cout << "error: " << (U.leftCols(r) * V.leftCols(r).transpose() - M).norm() / M.norm() << std::endl;
	
	int rowIndices[50];
	std::iota(rowIndices, rowIndices+50, 0);
	int rowLimits[] = {0, 50};
	int nRowClusters = 1;
	int colIndices[25];
	std::iota(colIndices, colIndices+25, 0);
	int colLimits[] = {0, 25};
	int nColClusters = 1;
	int nBlocks = 1;
	int rowBlocks[] = {0};
	int colBlocks[] = {0};
	
	double input[25], output[50];
	std::fill(input, input+25, 1.);
	
	std::fill(output, output+5, 0.);
	
	aca_multiply(M,
		rowIndices, rowLimits, nRowClusters,
		colIndices, colLimits, nColClusters,
		nBlocks, rowBlocks, colBlocks,
		input, output,
		1e-3, 50);
		
	std::ostream_iterator<double> out_it(std::cout, " ");
	std::copy(output, output+50, out_it);
	
	std::cout << "\n";
	
	// direct evaluation
	std::fill(output, output+50, 0.);
	for (int i = 0; i < 50; ++i)
		for (int j = 0; j < 25; ++j)
			output[rowIndices[i]] += M(rowIndices[i], colIndices[j]) * input[colIndices[j]];
	std::copy(output, output+50, out_it);

	return 0;
}


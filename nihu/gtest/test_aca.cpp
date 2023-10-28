#include "nihu/aca/aca.hpp"
#include <gtest/gtest.h>


// Assemble a random full-rank matrix and decompose into low-rank form
template <class T>
double ACA_Random_Tester(int r, int c, double eps, int maxRank, int &rank)
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(r, c), U(r, maxRank), V(c, maxRank);
	M.setRandom();
	rank = ACA::low_rank_approx(M, M.rows(), M.cols(), eps, maxRank, U, V);
	return (U.leftCols(rank) * V.leftCols(rank).transpose() - M).norm() / M.norm();
}

// Assemble a random low-rank matrix and decompose into low-rank form
template <class T>
double ACA_LowRank_Tester(int r, int c, int R, double eps, int &rank)
{
	int maxRank = std::min(r, c);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(r, c), U(r, maxRank), V(c, maxRank), u(r,1), v(c,1);
	M.setZero();
	for (int k = 0; k < R; ++k) {
		u.setRandom();
		v.setRandom();
		M += u * v.transpose();
	}
	rank = ACA::low_rank_approx(M, M.rows(), M.cols(), eps, maxRank, U, V);
	return (U.leftCols(rank) * V.leftCols(rank).transpose() - M).norm() / M.norm();
}


TEST(ACA, RandomReal) {
	int rank;
	double err = ACA_Random_Tester<double>(100, 300, 1e-5, 100, rank);
	EXPECT_NEAR(err, 0., 1e-3);
}

TEST(ACA, RandomComplex) {
	int rank;
	double err = ACA_Random_Tester<std::complex<double> >(100,300,1e-5,100, rank);
	EXPECT_NEAR(err, 0., 1e-3);
}


TEST(ACA, LowRankReal)
{
	int rank;
	for (int R = 1; R < 100; ++R)
	{
		double err = ACA_LowRank_Tester<double>(100, 300, R, 1e-10, rank);
		EXPECT_LE(rank, R);
		EXPECT_NEAR(err, 0., 1e-10);
	}
}

TEST(ACA, LowRankComplex)
{
	int rank;
	for (int R = 1; R < 100; ++R)
	{
		double err = ACA_LowRank_Tester<std::complex<double> >(100, 300, R, 1e-10, rank);
		EXPECT_LE(rank, R);
		EXPECT_NEAR(err, 0., 1e-10);
	}
}


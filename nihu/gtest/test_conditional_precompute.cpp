#include <Eigen/Core>
#include "util/conditional_precompute.hpp"
#include <gtest/gtest.h>

typedef Eigen::Matrix<double, 20, 3> mtx_type;

class Mtx
{
	static unsigned num_evals;

public:
	static mtx_type eval(void)
	{
		num_evals++;
		return mtx_type::Random();
	}
	
	static unsigned get_num_evals(void) { return num_evals; }
	static void reset(void) { num_evals = 0; }
};

unsigned Mtx::num_evals = 0;

typedef NiHu::conditional_precompute<NiHu::matrix_function_complexity::general, Mtx> general_type;
static_assert(std::is_same<mtx_type, general_type::return_type>::value, "general conditional_precompute should return the matrix value type");

typedef NiHu::conditional_precompute<NiHu::matrix_function_complexity::constant, Mtx> constant_type;
static_assert(std::is_same<mtx_type const &, constant_type::return_type>::value, "constant conditional_precompute should return the matrix constant reference");

typedef NiHu::conditional_precompute<NiHu::matrix_function_complexity::zero, Mtx> zero_type;
static_assert(std::is_same<decltype(mtx_type::Zero()), zero_type::return_type>::value, "zero conditional_precompute should return the zero matrix type");

TEST(ConditionalPrecompute, Static)
{
	Mtx::reset();
	auto g1 = general_type::eval();
	auto g2 = general_type::eval();
	EXPECT_EQ(Mtx::get_num_evals(), 2);
	EXPECT_NE(g1, g2);
	
	Mtx::reset();
	auto c1 = constant_type::eval();
	auto c2 = constant_type::eval();
	EXPECT_EQ(Mtx::get_num_evals(), 1);
	EXPECT_EQ(c1, c2);
	
	Mtx::reset();
	auto z1 = zero_type::eval();
	auto z2 = zero_type::eval();
	EXPECT_EQ(Mtx::get_num_evals(), 0);
	EXPECT_EQ(z1, z2);
	EXPECT_EQ(z1, mtx_type::Zero());
}


typedef NiHu::conditional_precompute_instance<NiHu::matrix_function_complexity::general, Mtx> general_type_instance;
typedef NiHu::conditional_precompute_instance<NiHu::matrix_function_complexity::constant, Mtx> constant_type_instance;
typedef NiHu::conditional_precompute_instance<NiHu::matrix_function_complexity::zero, Mtx> zero_type_instance;

TEST(ConditionalPrecompute, Instance)
{
	Mtx::reset();
	general_type_instance g;
	auto g1 = g();
	auto g2 = g();
	EXPECT_EQ(Mtx::get_num_evals(), 2);
	EXPECT_NE(g1, g2);
	
	Mtx::reset();
	constant_type_instance c;
	auto c1 = c();
	auto c2 = c();
	EXPECT_EQ(Mtx::get_num_evals(), 1);
	EXPECT_EQ(c1, c2);
	
	Mtx::reset();
	zero_type_instance z;
	auto z1 = z();
	auto z2 = z();
	EXPECT_EQ(Mtx::get_num_evals(), 0);
	EXPECT_EQ(z1, z2);
	EXPECT_EQ(z1, mtx_type::Zero());
}



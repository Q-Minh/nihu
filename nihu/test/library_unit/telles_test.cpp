/**
 * @file telles_test.cpp
 * @brief Tests for the Telles transformation
 * @ingroup library_test
 * 
 * @details
 * Some examples from:
 * L. Jun, G. Beer, J.L. Meek: Efficient evaluation of integrals of order 
 *   1/r 1/r^2, 1/r^3 using Gauss quadrature. Engineering Analysis Vol. 2,
 *   Issue 3, pp. 118-128 (1985) DOI: 10.1016/0264-682X(85)90014-0
 */

#include "core/gaussian_quadrature.hpp"
#include "library/telles_1987.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

typedef NiHu::gaussian_quadrature<NiHu::line_domain> quad_1d_t;
typedef NiHu::quadrature_traits<quad_1d_t>::domain_t::xi_t xi_1d_t;

typedef NiHu::gaussian_quadrature<NiHu::quad_domain> quad_2d_t;
typedef NiHu::quadrature_traits<quad_2d_t>::domain_t::xi_t xi_2d_t;


template <class QuadDerived, class func>
double integrate(func f, NiHu::quadrature_base<QuadDerived> const& q)
{
	double result = 0;
	for (auto it = q.begin(); it != q.end(); ++it) {
		result += it->get_w() * f(it->get_xi());
	}
	return result;
}

void show_results(double I_ana, double I_gauss, double I_telles, size_t n_gauss, size_t n_telles)
{
	std::cout << "\tAnalytical:\t\t" << std::setprecision(12) << std::fixed << I_ana << std::endl;
	std::cout << "\tGauss (" << n_gauss << " points):\t" << std::setprecision(12) << std::fixed 
	          << I_gauss << " (rel. err.: " << std::setprecision(3) << std::scientific << (I_gauss - I_ana)/I_ana << ")" << std::endl;
	std::cout << "\tTelles (" << n_telles << " points):\t" << std::setprecision(12) << std::fixed 
	          << I_telles << " (rel. err.: " << std::setprecision(3) << std::scientific << (I_telles - I_ana)/I_ana << ")" << std::endl;
	std::cout << std::endl;
}

/* Test 1 */
double f1(xi_1d_t const& xi)
{
	return std::log(std::abs(0.3 + xi(0)));
}

void test1()
{
	quad_1d_t q(19);
	quad_1d_t q_t = NiHu::telles_transform(q, xi_1d_t(-0.3));
	
	std::cout << "Test 1: int(log|0.3 + x|)dx, x = -1 .. 1" << std::endl;
	double I_ana = -1.908598917;
	double I_gauss = integrate(f1, q);
	double I_telles = integrate(f1, q_t);
	
	show_results(I_ana, I_gauss, I_telles, q.size(), q_t.size());
}

/* Test 2 */
double f2(xi_1d_t const& xi)
{
	return 1.0 / ((1.1 - xi(0))*(1.1 - xi(0)));
}

void test2()
{
	quad_1d_t q(19);
	quad_1d_t q_t = NiHu::telles_transform(q, xi_1d_t(1.1));
	
	std::cout << "Test 2: int(1/(1.1-x)^2)dx, x = -1 .. 1" << std::endl;
	double I_ana = 9.52380952;
	double I_gauss = integrate(f2, q);
	double I_telles = integrate(f2, q_t);
	
	show_results(I_ana, I_gauss, I_telles, q.size(), q_t.size());
}

/* Test 3 */
double f3(xi_1d_t const& xi)
{
	return 1.0 / ((1.004 - xi(0))*(1.004 - xi(0)));
}

void test3()
{
	quad_1d_t q(19);
	quad_1d_t q_t = NiHu::telles_transform(q, xi_1d_t(1.004));
	
	std::cout << "Test 3: int(1/(1.004-x)^2)dx, x = -1 .. 1" << std::endl;
	double I_ana = 249.500998;
	double I_gauss = integrate(f3, q);
	double I_telles = integrate(f3, q_t);
	
	show_results(I_ana, I_gauss, I_telles, q.size(), q_t.size());
}

/* Test 4 */
double f4(xi_2d_t const& xi)
{
	return 1.0 / std::sqrt((1.004 - xi(0))*(1.004 - xi(0)) + (1.004 - xi(1))*(1.004 - xi(1)));
}

void test4()
{
	quad_2d_t q(15);
	quad_2d_t q_t = NiHu::telles_transform(q, xi_2d_t(1.004, 1.004));
	
	std::cout << "Test 4: int(1/sqrt((1.004-x)^2+(1.004-y)^2))dxdy, x,y = -1 .. 1" << std::endl;
	double I_ana = 3.476318;
	double I_gauss = integrate(f4, q);
	double I_telles = integrate(f4, q_t);
	
	show_results(I_ana, I_gauss, I_telles, q.size(), q_t.size());
}

int main()
{
	test1();
	test2();
	test3();
	test4();
	return 0;
}
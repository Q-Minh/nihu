// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/math/constants/constants.hpp>

#include "nihu/core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
#include "nihu/library/helmholtz_kernel.hpp"
#include "nihu/library/helmholtz_singular_integrals.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

/// \todo create circle again
NiHu::mesh<tmp::vector<NiHu::line_1_elem> > create_circle(double R, int N)
{
	using namespace boost::math::double_constants;
	
	dMatrix nodes(N, 2);
	uMatrix elements(N,3);
	double dphi = two_pi / N;
	for (int i = 0; i < N; ++i)
	{
		double phi = i * dphi;
		nodes(i,0) = R * std::cos(phi);
		nodes(i,1) = R * std::sin(phi);
		elements(i,0) = NiHu::line_1_elem::id;
		elements(i,1) = i % N;
		elements(i,2) = (i+1) % N;
	}
	return NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
}

int main(void)
{
	double r = 1;
	int N = 200;
	auto mesh = create_circle(r, N);

	auto const &trial_space = NiHu::constant_view(mesh);
	auto const &test_space = trial_space;
	
	cMatrix L(N, N), M(N, N), I(N,N);
	L.setZero(); M.setZero(); I.setZero();
	
	double k = 1.0;

	std::complex<double> q0(1.0);
	std::complex<double> jkr(0,k*r);
	std::complex<double> panal(-r / (1.+jkr));
	
	auto I_op = NiHu::identity_integral_operator();
	auto L_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_SLP_kernel<double>(k));
	auto M_op = NiHu::create_integral_operator(NiHu::helmholtz_2d_DLP_kernel<double>(k));

	I << test_space * I_op[trial_space];
	L << test_space * L_op[trial_space];
	M << test_space * M_op[trial_space];
	
	cMatrix q(N,1);
	q.setConstant(q0);
	
	cMatrix p = (M - .5 * I).colPivHouseholderQr().solve(L * q);
	
	std::cout << "Analytic: " << panal << std::endl;
	std::cout << "Mean numeric: " << p.mean() << std::endl;
	std::cout << "log10 Error: " << std::log10((p.array()-panal).matrix().norm() / std::abs(panal) / std::sqrt(N)) << std::endl;

	return 0;
}

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

#include "core/weighted_residual.hpp"
#include "util/mex_matrix.hpp"
#include "library/elastodynamics_kernel.hpp"
#include "library/elastodynamics_singular_integrals.hpp"
#include "library/lib_element.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;
typedef NiHu::mex::complex_matrix<double> cMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	using namespace boost::math::double_constants;
	
	dMatrix surf_nodes(rhs[0]), surf_elem(rhs[1]);
	auto surf_mesh = NiHu::create_mesh(surf_nodes, surf_elem, NiHu::quad_1_tag());
	auto const &surf_sp = NiHu::constant_view(surf_mesh, NiHu::_3d());
	size_t n = surf_sp.get_num_dofs();

	dMatrix field_nodes(rhs[2]), field_elem(rhs[3]);
	auto field_mesh = NiHu::create_mesh(field_nodes, field_elem, NiHu::quad_1_tag());
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field_mesh, NiHu::_3d()));
	size_t m = field_sp.get_num_dofs();

	cMatrix Us(n, n, lhs[0]), Ts(n, n, lhs[1]), Uf(m, n, lhs[2]), Tf(m, n, lhs[3]);
	
	double nu = NiHu::mex::get_scalar<double>(rhs[4]);
	double rho = NiHu::mex::get_scalar<double>(rhs[5]);
	double mu = NiHu::mex::get_scalar<double>(rhs[6]);
	double freq = NiHu::mex::get_scalar<double>(rhs[7]);

	auto U_op = NiHu::create_integral_operator(NiHu::elastodynamics_3d_U_kernel(nu, rho, mu, two_pi * freq));
	auto T_op = NiHu::create_integral_operator(NiHu::elastodynamics_3d_T_kernel(nu, rho, mu, two_pi * freq));

	Us << dirac(surf_sp) * U_op[surf_sp];
	Ts << dirac(surf_sp) * T_op[surf_sp];

	Uf << field_sp * U_op[surf_sp];
	Tf << field_sp * T_op[surf_sp];
}


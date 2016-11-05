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

#include "core/weighted_residual.hpp"
#include "library/covariance_kernel.hpp"
#include "util/mex_matrix.hpp"
#include "library/lib_element.hpp"
#include <Eigen/Eigenvalues>

typedef mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	dMatrix nodes(rhs[0]), elem(rhs[1]);
	auto mesh = create_mesh(nodes, elem, line_1_volume_tag());

	auto const &w = constant_view(mesh);

	double sigma = *mxGetPr(rhs[2]);
	double d = *mxGetPr(rhs[3]);
	auto C = create_integral_operator(covariance_kernel(sigma, d));
	auto I = identity_integral_operator();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> D(w.get_num_dofs(), w.get_num_dofs()), B(w.get_num_dofs(), w.get_num_dofs());
	D.setZero();
	B.setZero();
	
	D << (w * C[w]);
	B << (w * I[w]);
	
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > es(D, B);

	dMatrix G(D.rows(), D.rows(), lhs[0]);
	dMatrix lambda(D.rows(), 1, lhs[1]);
	
	for (int i = 0; i < D.rows(); ++i)
	{
		lambda(i) = es.eigenvalues()(i);
		for (int j = 0; j < D.rows(); ++j)
			G(i,j) = es.eigenvectors()(i,j);
	}
}


#include <fstream>
#include <iostream>
#include <string>

#include "interface/read_off_mesh.hpp"

#include "library/lib_element.hpp"
#include "library/elastostatics"

#include "core/function_space.hpp"
#include "core/integral_operator.hpp"
#include "core/weighted_residual.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;

int main(int argc, char const* argv[])
{
	try
	{
		// call: prog.exe meshname xctname respname
		std::string meshname(argv[1]);
		std::string excname(argv[2]);
		std::string respname(argv[3]);

		auto mesh = NiHu::read_off_mesh(meshname, NiHu::quad_1_tag());
		std::cout << "Mesh read. Number of elements: " << mesh.get_num_elements() << std::endl;

		// create function spaces
		auto const& trial_space = NiHu::constant_view(mesh, NiHu::_3d());
		auto const& test_space = NiHu::dirac(trial_space);
		size_t n = trial_space.get_num_dofs();
		size_t m = test_space.get_num_dofs();
		std::cout << "Function spaces created. Number of DOFs: " << m << " x " << n << std::endl;

		// create kernel
		double const nu = .33;
		double const mu = 1e9;
		NiHu::elastostatics_3d_U_kernel ukernel(nu, mu);
		NiHu::elastostatics_3d_T_kernel tkernel(nu, mu);
		auto U_op = NiHu::create_integral_operator(ukernel);
		auto T_op = NiHu::create_integral_operator(tkernel);


		// read excitation
		dVector t(n, 1), rhs, u;

		{
			std::ifstream ifs(excname);
			if (!ifs)
				throw std::runtime_error("Could not open file " + excname);
			for (size_t i = 0; i < n; ++i)
				ifs >> t(i, 0);
		}

		{
			dMatrix U(n, n);
			U.setZero();
			U << test_space * U_op[trial_space];
			rhs = U * t;
		}

		{
			dMatrix T = -.5 * dMatrix::Identity(n, n);
			T << test_space * T_op[trial_space];
			u = T.colPivHouseholderQr().solve(rhs);
		}

		{
			std::ofstream ofs(respname);
			if (!ofs)
				throw std::runtime_error("Could not open file " + respname);
			for (size_t i = 0; i < n; ++i)
				ofs << u(i, 0) << std::endl;
		}
	}
	catch (std::exception & exc)
	{
		std::cerr << "Standard exception caught: " << exc.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unhandled exception caught" << std::endl;
	}

	return 0;
}

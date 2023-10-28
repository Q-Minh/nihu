#include <iostream>
#include <fstream>
#include <sstream>

#include <omp.h>

#include "core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/quad_1_gauss_field.hpp"

#include<Eigen/IterativeLinearSolvers>

#define NUM_PROCESSORS 8

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cVector;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

void read_off_data(std::string const &fname, dMatrix &nodes, uMatrix &elements)
{
	// open mesh file for reading
	std::ifstream is(fname);
	if (!is)
		throw std::runtime_error("Error reading mesh file");

	// read header from file (first row is 'OFF')
	std::string header;
	if (!(is >> header) || header != "OFF")
		throw std::runtime_error("Possibly invalid off file");

	// read number of nodes and number of elements, nEdges is dropped
	unsigned nNodes, nElements, nEdges;
	if (!(is >> nNodes >> nElements >> nEdges))
		throw std::runtime_error("Error reading number of mesh entries");

	// read nodes
	nodes.resize(nNodes, 3);
	for (unsigned i = 0; i < nNodes; ++i)
		if (!(is >> nodes(i,0) >> nodes(i,1) >> nodes(i,2)))
			throw std::runtime_error("Error reading mesh nodes");

	// read elements
	elements.resize(nElements, 5);
	for (unsigned i = 0; i < nElements; ++i)
	{
		unsigned nvert;
		if (!(is >> nvert))
			throw std::runtime_error("Error reading mesh elements");
		for (unsigned c = 0; c < nvert; ++c)
			if (!(is >> elements(i,c+1)))
				throw std::runtime_error("Error reading mesh elements");
		elements(i,0) = NiHu::quad_1_elem::id;
	}

	is.close();
}




template <class TestSpace, class TrialSpace, class FieldSpace>
void solve(TestSpace const &test, TrialSpace const &trial, FieldSpace const &field_sp,
	char const *pattern)
{
	size_t nDof = trial.get_num_dofs();
	size_t M = field_sp.get_num_dofs();
	
	double const c = 340.;
	double const rho = 1.3;
	std::complex<double> const J(0., 1.);
	
	size_t nFreqs = 120;
	size_t nBlock = nFreqs / NUM_PROCESSORS;
	
	dMatrix fvec(nBlock, NUM_PROCESSORS);
	for (size_t i = 0; i < nBlock; ++i)
		for (size_t j = 0; j < NUM_PROCESSORS; ++j)
			fvec(i,j) = (NUM_PROCESSORS*i+j) * .5;
	
#pragma omp parallel for num_threads(NUM_PROCESSORS)
	for (size_t i = 0; i < nFreqs; ++i)
	{
		double f = fvec(i);
		double om = 2.*M_PI*f;
		double k = om / c;
		
#ifdef BM
		std::complex<double> alpha(0., 1./k);
#endif		
		
		// create integral operators
		auto Gop = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
		auto Hop = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
		auto Iop = NiHu::identity_integral_operator();
#ifdef BM
		auto Htop = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
		auto Dop = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));
#endif
		
		// create excitation
		cVector qs(nDof, 1);
		double v0 = 1e-3;
		qs.setConstant(-J*om*rho*v0);
		
		// compute rhs
		cVector rhs;
		
		{
			cMatrix Gs(nDof, nDof);
			Gs.setZero();
			std::cout << "Integrating Gs at f = " << f << std::endl;
			Gs << test * Gop[trial];
			
#ifdef BM
			std::cout << "Integrating Hts at f = " << f << std::endl;
			Gs << test * (alpha * Htop)[trial];
			Gs << test * (alpha/2. * Iop)[trial];
#endif			
			
			rhs = Gs * qs;
		}
		
		cVector ps;
		size_t num_iters = 0;
		
		{
			// create matrices
			cMatrix Hs(nDof, nDof);
			Hs.setZero();
			
#ifdef BM
			std::cout << "Integrating Ds at f = " << f << std::endl;
			Hs << test * (alpha * Dop)[trial];
#endif

			std::cout << "Integrating Hs at f = " << f << std::endl;
			Hs << test * (-.5 * Iop)[trial];
			Hs << test * Hop[trial];
			
			// solve linear system
			std::cout << "Solving linear system" << std::endl;
#ifdef ITERATIVE
			Eigen::BiCGSTAB<cMatrix> solver(Hs);
			solver.setTolerance(1e-8);
			ps = solver.solve(rhs);
			num_iters = solver.iterations();
			
			std::cout << "f:               " << f << '\n';
			std::cout << "#iterations:     " << solver.iterations() << '\n';
			std::cout << "estimated error: " << solver.error()      << std::endl;
#else
			ps = Hs.colPivHouseholderQr().solve(rhs);
			std::cout << "f:               " << f << '\n';
#endif
		}
		
		{
			// export ps
			std::stringstream ss;
			ss << pattern << "_" << f << "ps.res";
			std::ofstream os(ss.str().c_str());
			os << f << ' ' << num_iters << std::endl;
			for (size_t i = 0; i < nDof; ++i)
				os << ps(i).real() << ' ' << ps(i).imag() << std::endl;
			os.close();
		}
		
		cVector pf;
		
		{
			cMatrix Gf(M, nDof), Hf(M, nDof);
			Gf.setZero();
			Hf.setZero();

			std::cout << "Integrating Gf" << std::endl;
			Gf << field_sp * Gop[trial];
			std::cout << "Integrating Hf" << std::endl;
			Hf << field_sp * Hop[trial];
			
			pf = Hf * ps - Gf * qs;
		}
		
		{
			// export pf
			std::stringstream ss;
			ss << pattern << "_" << f << "pf.res";
			std::ofstream os(ss.str().c_str());
			os << f << std::endl;
			for (size_t i = 0; i < M; ++i)
				os << pf(i).real() << ' ' << pf(i).imag() << std::endl;
			os.close();
		}
	}
}

int main(int argc, char **argv)
{
	if (argc < 4)
		std::cerr << "usage: prog meshname fieldname pattern" << std::endl;
	
	char const *meshname = argv[1];
	char const *fieldname = argv[2];
	char const *pattern = argv[3];
	
#ifdef GAUSS
	// read mesh file
	uMatrix elements;
	dMatrix nodes;
	read_off_data(meshname, nodes, elements);

	// assemble field matrix
	size_t nElements = elements.rows();
	uMatrix fields(nElements, 1+4+4);
	for (size_t e = 0; e < nElements; ++e)
	{
		fields(e,0) = NiHu::quad_1_gauss_field::id;
		// assign nodes
		for (size_t c = 0; c < 4; ++c)
			fields(e,c+1) = elements(e,c+1);
		// assign DOF
		for (size_t c = 0; c < 4; ++c)
			fields(e,c+1+4) = 4*e+c;
	}
	
	// create function space
	auto trial = NiHu::create_function_space(nodes, fields, NiHu::quad_1_gauss_field_tag());
#else
	auto mesh = NiHu::read_off_mesh(meshname, NiHu::quad_1_tag(), NiHu::tria_1_tag());
	auto const &trial = NiHu::constant_view(mesh);
#endif
	auto const &test = NiHu::dirac(trial);
	
	// read mesh files
	auto field = NiHu::read_off_mesh(fieldname, NiHu::quad_1_tag(), NiHu::tria_1_tag());
	
	// create function space
	auto const &field_sp = NiHu::dirac(NiHu::constant_view(field));
	
	solve(test, trial, field_sp, pattern);
	
	return 0;
}

#include "../core/kernel.hpp"

template <class Space, class Layer>
class perturbed_laplace_kernel;


template <class Scalar>
class perturbed_laplace_kernel<space_3d<Scalar>, potential::SLP>
	: public kernel_base<perturbed_laplace_kernel<space_3d<Scalar>, potential::SLP> >
{
public:
	typedef kernel_base<perturbed_laplace_kernel<space_3d<Scalar>, potential::SLP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	result_t operator()(
		typename base_t::x_t const &x,
		typename base_t::y_t const &y,
		typename base_t::x_t const &ux,
		typename base_t::x_t const &uy,
		Scalar relative_jacobian) const
	{
		typename base_t::x_t rvec = y - x;
		auto r = rvec.norm();
		typename base_t::x_t evec = rvec.normalized();
		auto G = 1. / r / (4. * M_PI);
		auto gradG = -evec / (4. * M_PI) / r / r;
		auto dGtau = gradG.dot(uy - ux);
		return dGtau + G * relative_jacobian;
	}

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		return (*this)(test_input.get_x(), trial_input.get_x(),
			test_input.get_u(), trial_input.get_u(),
			trial_input.get_relative_jacobian());
	}
};

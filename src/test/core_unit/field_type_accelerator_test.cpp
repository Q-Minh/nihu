#include "bem/field_type_accelerator.hpp"
#include "bem/gaussian_quadrature.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

template <class Family, class Field>
struct tester
{
	struct type
	{
		typedef typename Field::domain_t domain_t;

		typedef quadrature_elem<
			typename domain_t::xi_t, typename domain_t::scalar_t
		> qe_t;

		typedef field_type_accelerator_elem<
			Field, Family, acceleration::soft
		> soft_accel_t;

		typedef field_type_accelerator_elem<
			Field, Family, acceleration::hard
		> hard_accel_t;

		void operator() (void)
		{
			// instantiate a quadrature
			typename quadrature_type<Family, domain_t>::type q(5);

			std::cout << "soft accelerator:\n";

			auto a_soft = static_cast<soft_accel_t const &>(q[0]);
			std::cout << a_soft.get_xi().transpose() << '\t' <<
				a_soft.get_w() << '\t' <<
				a_soft.get_N().transpose() << '\n';

			std::cout << "hard accelerator:\n";

			hard_accel_t a_hard(q[0]);
			std::cout << a_hard.get_xi().transpose() << '\t' <<
				a_hard.get_w() << '\t' <<
				a_hard.get_N().transpose() << '\n';

			field_type_accelerator<Field, Family, acceleration::hard> acc_hard(q);
			for (auto it = acc_hard.begin(); it != acc_hard.end(); ++it)
				std::cout << it->get_w() << std::endl;

			auto acc_soft = static_cast<field_type_accelerator<Field, Family, acceleration::soft> >(q);
			for (auto it = acc_soft.begin(); it != acc_soft.end(); ++it)
				std::cout << it->get_w() << std::endl;
		}
	};
};


int main(void)
{
	std::cout << "Testing regular quadratures" << std::endl << std::endl;
	tmp::call_each<
		tmp::vector<field_view<quad_1_elem, field_option::isoparametric> >,
		tester<gauss_family_tag, tmp::_1>
	>();

	return 0;
}


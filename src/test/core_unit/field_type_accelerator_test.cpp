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
			field_type_accelerator<Field, Family, acceleration::hard> ha(5);

			// instantiate a quadrature
			typename quadrature_type<Family, domain_t>::type q(5);

			std::cout << "soft accelerator:\n";

			auto const &a_soft = static_cast<soft_accel_t const &>(q[0]);
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

			std::cout << std::endl;

			auto const &acc_soft = static_cast<field_type_accelerator<Field, Family, acceleration::soft> const &>(q);
			for (auto it = acc_soft.begin(); it != acc_soft.end(); ++it)
				std::cout << it->get_w() << std::endl;
		}
	};
};



template <class Family, class Field>
struct dirac_tester
{
	struct type
	{
		void operator() (void)
		{
			field_type_accelerator<Field, Family, acceleration::hard> dirac_acc;

			for (auto it = dirac_acc.begin(); it != dirac_acc.end(); ++it)
				std::cout << it->get_w() << '\t' <<
					it->get_N().transpose() << std::endl;

			std::cout << std::endl;
		}
	};
};


template <class TestField, class TrialField>
struct test_dual_regular
{
	static unsigned const MAX_ORDER = 10;

	typedef store<field_type_accelerator_pool<
		TestField, gauss_family_tag, acceleration::soft, MAX_ORDER
	> > test_store_t;

	typedef store<field_type_accelerator_pool<
		TrialField, gauss_family_tag, acceleration::soft, MAX_ORDER
	> > trial_store_t;

	struct type
	{
		void operator()(void)
		{
			auto dual_acc = create_dual_field_type_accelerator(
				test_store_t::m_data[4], trial_store_t::m_data[4], iteration::diadic());
			for (auto it = dual_acc.begin(); it != dual_acc.end(); ++it)
			{
				std::cout << it.get_first()->get_w() * it.get_second()->get_w() << std::endl;
			}
		}
	};
};





int main(void)
{
	tmp::d_call_each<
		tmp::vector<field_view<quad_1_elem, field_option::constant> >,
		tmp::vector<field_view<quad_1_elem, field_option::constant> >,
		test_dual_regular<tmp::_1, tmp::_2>
	>();

	return 0;
}


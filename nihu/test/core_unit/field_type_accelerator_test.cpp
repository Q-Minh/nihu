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

#include "nihu/core/field_type_accelerator.hpp"
#include "nihu/core/gaussian_quadrature.hpp"
#include "nihu/util/store_pattern.hpp"
#include "nihu/tmp/control.hpp"
#include "nihu/tmp/vector.hpp"
#include "nihu/library/lib_element.hpp"

template <class Family, class Field>
struct tester
{
	struct type
	{
		typedef typename Field::domain_t domain_t;

		typedef NiHu::quadrature_elem<
			typename domain_t::xi_t, typename domain_t::scalar_t
		> qe_t;

		typedef NiHu::field_type_accelerator_elem<
			Field, NiHu::gauss_family_tag, NiHu::acceleration::soft
		> soft_accel_t;

		typedef NiHu::field_type_accelerator_elem<
			Field, NiHu::gauss_family_tag, NiHu::acceleration::hard
		> hard_accel_t;

		void operator() (void)
		{
			NiHu::field_type_accelerator<Field, NiHu::gauss_family_tag, NiHu::acceleration::hard> ha(5);

			// instantiate a quadrature
			typename NiHu::quadrature_type<NiHu::gauss_family_tag, domain_t>::type q(5);

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

			NiHu::field_type_accelerator<Field, NiHu::gauss_family_tag, NiHu::acceleration::hard> acc_hard(q);
			for (auto it = acc_hard.begin(); it != acc_hard.end(); ++it)
				std::cout << it->get_w() << std::endl;

			std::cout << std::endl;

			auto const &acc_soft = static_cast<NiHu::field_type_accelerator<Field, NiHu::gauss_family_tag, NiHu::acceleration::soft> const &>(q);
			for (auto it = acc_soft.begin(); it != acc_soft.end(); ++it)
				std::cout << it->get_w() << std::endl;
		}
	};
};



template <class Field>
struct test_dirac
{
	typedef NiHu::dirac_field<Field> dfield_t;
	struct type
	{
		void operator() (void)
		{
			// test field type accelerator
			NiHu::field_type_accelerator<dfield_t, NiHu::gauss_family_tag, NiHu::acceleration::soft> dirac_acc(0);

			for (auto it = dirac_acc.begin(); it != dirac_acc.end(); ++it)
				std::cout << it->get_w() << '\t'
					<< it->get_N().transpose() << std::endl;

			std::cout << std::endl;
		}
	};
};


template <class TestField, class TrialField>
struct test_dual_regular
{
	static unsigned const MAX_ORDER = 10;

	typedef NiHu::store<NiHu::field_type_accelerator_pool<
		TestField, NiHu::gauss_family_tag, NiHu::acceleration::soft, MAX_ORDER
	> > test_store_t;

	typedef NiHu::store<NiHu::field_type_accelerator_pool<
		TrialField, NiHu::gauss_family_tag, NiHu::acceleration::soft, MAX_ORDER
	> > trial_store_t;

	struct type
	{
		void operator()(void)
		{
			auto dual_acc = NiHu::create_dual_field_type_accelerator(
				test_store_t::get_data()[4], trial_store_t::get_data()[4], NiHu::iteration::diadic());
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
		tmp::vector<NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> >,
		tmp::vector<NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> >,
		test_dual_regular<tmp::_1, tmp::_2>
	>();

	tmp::call_each<
		tmp::vector<NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> >,
		test_dirac<tmp::_1>
	>();

	return 0;
}


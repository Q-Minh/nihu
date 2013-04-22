#include "../bem/gaussian_quadrature.hpp"
#include "../bem/duffy_quadrature.hpp"
#include "../tmp/control.hpp"

template <class Family, class Domain>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			typename quadrature_type<Family, Domain>::type q(5);
			std::cout << q << std::endl;
		}
	};
};


void test_transform(void)
{
	gauss_quad q(4);
	Eigen::Matrix<double, 4, 2> coords;
	coords <<
		0.0, 0.0,
		1.0, 0.0,
		1.0, 1.0,
		1.0, 1.0;
	std::cout << q.transform<quad_1_shape_set>(coords) << std::endl;
}


template <class Family, class LSet>
struct duffy_test
{
	struct type {
		void operator()(void)
		{
			typedef duffy_quadrature<Family, LSet> duffy_t;

			auto d_corner = duffy_t::on_corner(3, 2);
			std::cout << d_corner << std::endl;

			auto d_face = duffy_t::on_face(3, LSet::domain_t::get_center());
			std::cout << d_face << std::endl;
		}
	};
};


int main(void)
{
	std::cout << "Testing regular quadratures" << std::endl << std::endl;
	tmp::call_each<
		tmp::vector<line_domain, tria_domain, quad_domain>,
		tester<gauss_family_tag, tmp::_1>
	>();

	std::cout << "Testing quadrature transforms" << std::endl << std::endl;
	test_transform();

	std::cout << "Testing duffy quadratures" << std::endl << std::endl;
	tmp::call_each<
		tmp::vector<tria_1_shape_set, quad_1_shape_set>,
		duffy_test<gauss_family_tag, tmp::_1>
	>();

	return 0;
}


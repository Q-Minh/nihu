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
			std::cout << "Testing regular quadratures" << std::endl << std::endl;

			typename quadrature_type<Family, Domain>::type q(5);
			std::cout << q << std::endl;
		}
	};
};

void test_transform(void)
{
	std::cout << "Testing quadrature transforms" << std::endl << std::endl;

	gauss_quad q(4);
	Eigen::Matrix<double, 4, 2> coords;
	coords <<
		0.0, 0.0,
		1.0, 0.0,
		1.0, 1.0,
		1.0, 1.0;
	std::cout << q.transform<quad_1_shape_set>(coords) << std::endl;
}



void test_duffy(void)
{
	std::cout << "Testing duffy quadratures" << std::endl << std::endl;

	typedef quad_1_shape_set lset_t;
	typedef duffy_quadrature<gauss_family_tag, lset_t> duffy_t;

	auto d_corner = duffy_t::on_corner(3, 2);
	std::cout << d_corner << std::endl;

	auto d_face = duffy_t::on_face(3, duffy_t::xi_t(.0, .0));
	std::cout << d_face << std::endl;
}

int main(void)
{
	tmp::call_each<
		tmp::vector<line_domain, tria_domain, quad_domain>,
		tester<gauss_family_tag, tmp::_1>
	>();
	test_transform();
	test_duffy();

	return 0;
}


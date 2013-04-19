#include "../bem/quadrature.hpp"

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
	std::cout << gauss_quad::singular_quadrature_inside(5, quad_domain::get_center());
	std::cout << gauss_tria::singular_quadrature_inside(5, tria_domain::get_center());
}


#include "../bem/duffy_quadrature.hpp"

void test_duffy(void)
{
	typedef tria_1_shape_set lset_t;
	typedef quadrature_type<gauss_family_tag, typename lset_t::domain_t>::type	duffy_t;
	duffy_t duf = duffy_quadrature_corner<lset_t>(5, 0);
	std::cout << duf << std::endl;
}


int main(void)
{
	typedef gauss_family_tag gauss;
	typedef tmp::vector<line_domain, tria_domain, quad_domain> DomainSequence;
	//typedef tmp::vector<gauss_quad, gauss_tria> QuadSequence;

	tmp::call_each<DomainSequence, tester<gauss, tmp::_1> >();

	test_transform();

	test_duffy();

	return 0;
}


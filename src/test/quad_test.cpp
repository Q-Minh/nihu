#include "../bem/quadrature.hpp"

#include "../tmp/control.hpp"

template <class Family, class Domain>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			typename quadrature_domain_traits<Family, Domain>::quadrature_type q(5);
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


int main(void)
{
	typedef gauss_family_tag gauss;
	typedef tmp::vector<line_domain, tria_domain, quad_domain> DomainSequence;
	//typedef tmp::vector<gauss_quad, gauss_tria> QuadSequence;

	tmp::call_each<DomainSequence, tester<gauss, tmp::_1> >();

	test_transform();

	return 0;
}


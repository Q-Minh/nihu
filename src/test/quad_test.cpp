#include "../bem/gaussian_quadrature.hpp"

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


#include "../bem/duffy_quadrature.hpp"

void test_duffy(void)
{
	std::cout << "Testing duffy quadratures" << std::endl << std::endl;

	typedef tria_1_shape_set lset_t;
	typedef quadrature_type<gauss_family_tag, typename lset_t::domain_t>::type	duffy_t;

	duffy_t duf_corner = duffy_quadrature_corner<gauss_family_tag, lset_t>(3, 2);
	std::cout << duf_corner << std::endl;

	duffy_t duf_face = duffy_quadrature_face<gauss_family_tag, lset_t>(3, lset_t::xi_t(0.0, 0.0));
	std::cout << duf_face << std::endl;
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


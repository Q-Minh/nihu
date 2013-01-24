#include "quadrature.hpp"

#include "../tmp/control.hpp"

template <class D>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			gauss_quadrature<D> q(5);
			std::cout << q << std::endl;
			std::cout << "Sum of weights: " <<
				std::accumulate(q.begin(), q.end(), 0.0, [] (double x, quadrature_elem<D> &qe) {
				return x + qe.get_w();
			}) << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<line_domain, tria_domain, quad_domain> DomSequence;

	tmp::call_each<DomSequence, tester<tmp::_1> >();

	return 0;
}


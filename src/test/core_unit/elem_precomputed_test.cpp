#include "tmp/sequence.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include "core/element.hpp"
#include "library/lib_element.hpp"

#include <chrono>

template <class ElemType>
struct tester
{
	struct type
	{
		void operator()(void)
		{
			std::cout << "\n" << element_traits::name<ElemType>::value << " (" << ElemType::id << ")\n"
				<< "==========" << std::endl;
			const unsigned  N = 1e6;
			std::vector<ElemType> elements;
			elements.reserve(N);
			// construct elements

			typename ElemType::xi_t xi = ElemType::domain_t::get_center();

			auto tic = std::chrono::system_clock::now();
			for (unsigned i = 0; i < N; ++i)
				elements.push_back(ElemType(ElemType::coords_t::Random()));
			auto toc = std::chrono::system_clock::now();
			auto const_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();


			tic = std::chrono::system_clock::now();
			double n = 0.0;
			for (auto it = elements.begin(); it != elements.end(); ++it)
				n += it->get_normal(xi).norm();
			toc = std::chrono::system_clock::now();
			auto comp_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();

			std::cout << "construction time: " << const_elapsed << std::endl;
			std::cout << "computation time: " << comp_elapsed << std::endl;
			std::cout << "total time: " << const_elapsed + comp_elapsed << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			line_1_elem, line_2_elem,
			tria_1_elem, quad_1_elem,
			tria_2_elem, quad_2_elem, quad_28_elem
		>,
		tester<tmp::_1>
	>();
	return 0;
}

#include "../bem/shapeset.hpp"

#include "../tmp/control.hpp"

template <class shape_set>
struct shape_tester
{
	struct type {
		void operator()(void)
		{
			std::cout << "corners: " << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << it->transpose() << std::endl;
			std::cout << "shape values" << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << shape_set::eval_shape(*it).transpose() << std::endl;
			std::cout << "shape derivative values" << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << shape_set::eval_dshape(*it).transpose() << "sum: " << shape_set::eval_dshape(*it).sum() << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			tria_2_shape_set,
			quad_0_shape_set,
			quad_1_shape_set,
			quad_2_shape_set
		>,
		shape_tester<tmp::_1>
	>();

	return 0;
}

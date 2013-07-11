
#include "../tmp/control.hpp"
#include "../bem/shapeset.hpp"

template <class shape_set>
struct tester
{
	struct type {
		void operator()(void)
		{
			std::cout << std::endl << "shape set id: " << shape_set::id << std::endl;

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
			line_0_shape_set, tria_0_shape_set, quad_0_shape_set, brick_0_shape_set,
			line_1_shape_set, tria_1_shape_set, quad_1_shape_set, brick_1_shape_set,
			parallelogram_shape_set, tria_2_shape_set, quad_2_shape_set
		>,
		tester<tmp::_1>
	>();

	return 0;
}


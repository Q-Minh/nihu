#include "bem/interval.hpp"
#include "tmp/control.hpp"
#include <iostream>

typedef tmp::vector<
	break_point<std::ratio<2>, tmp::int_<6> >,
	break_point<ratio_infinite, tmp::int_<1> >
> inter1_t;

typedef tmp::vector<
	break_point<std::ratio<3>, tmp::int_<6> >,
	break_point<ratio_infinite, tmp::int_<2> >
> inter2_t;

template <class BP>
struct print_break_point
{
	struct type
	{
		void operator()(void)
		{
			std::cout << BP::x_value() << " - " << BP::y::value << std::endl;
		}
	};
};

int main(void)
{
	std::cout << "first interval: " << std::endl;
	tmp::call_each<inter1_t, print_break_point<tmp::_1> >();
	std::cout << "second interval: " << std::endl;
	tmp::call_each<inter2_t, print_break_point<tmp::_1> >();

	std::cout << "merged interval: " << std::endl;
	typedef merge_intervals<inter1_t, inter2_t>::type merged;
	tmp::call_each<merged, print_break_point<tmp::_1> >();

	std::cout << "some tests: " << std::endl;
	std::cout << .2 << ' ' << eval_interval<merged>(.2) << std::endl;
	std::cout << 2.2 << ' ' << eval_interval<merged>(2.2) << std::endl;
	std::cout << 4.2 << ' ' << eval_interval<merged>(4.2) << std::endl;
	return 0;
}

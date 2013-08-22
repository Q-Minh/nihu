#include "bem/interval.hpp"
#include "tmp/control.hpp"
#include <iostream>

typedef tmp::vector<
	break_point<std::ratio<0,1>, tmp::int_<6> >,
	break_point<std::ratio<2,1>, tmp::int_<3> >,
	break_point<std::ratio<5,2>, tmp::int_<1> >
> func_t;

template <class BP>
struct print_break_point
{
	struct type
	{
		void operator()(void)
		{
			std::cout << BP::x_value() << '-' << BP::y::value << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<func_t, print_break_point<tmp::_1> >();
	return 0;
}

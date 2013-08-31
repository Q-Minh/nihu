#include "core/interval.hpp"
#include "tmp/control.hpp"
#include <iostream>

typedef tmp::vector<
	break_point<std::ratio<2>, tmp::int_<8> >,
	break_point<ratio_infinite, tmp::int_<1> >
> inter1_t;

typedef tmp::vector<
	break_point<std::ratio<3>, tmp::int_<6> >,
	break_point<std::ratio<5>, tmp::int_<4> >,
	break_point<ratio_infinite, tmp::int_<2> >
> inter2_t;

typedef merge_intervals<inter1_t, inter2_t>::type merged_t;

template <class BP>
struct print_break_point
{
	struct type
	{
		void operator()(void)
		{
			std::cout << BP::x_value() << " : " << BP::y::value << std::endl;
		}
	};
};

int main(void)
{
	std::cout << "Test intervals" << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "I1 :" << std::endl;
	tmp::call_each<inter1_t, print_break_point<tmp::_1> >();
	std::cout << "I2 :" << std::endl;
	tmp::call_each<inter2_t, print_break_point<tmp::_1> >();
	std::cout << std::endl;

	std::cout << "Testing interval merge" << std::endl;
	std::cout << "======================" << std::endl;
	
	std::cout << "Merged : " << std::endl;
	tmp::call_each<merged_t, print_break_point<tmp::_1> >();
	std::cout << std::endl;

	std::cout << "Testing interval evaluation" << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << "  x : I1, I2, Merged " << std::endl;
	std::cout <<  .2 << " : " <<' ' << eval_interval<inter1_t>( .2) << ' ' << eval_interval<inter2_t>( .2) << ' ' << eval_interval<merged_t>( .2) << std::endl;
	std::cout << 2.2 << " : " <<' ' << eval_interval<inter1_t>(2.2) << ' ' << eval_interval<inter2_t>(2.2) << ' ' << eval_interval<merged_t>(2.2) << std::endl;
	std::cout << 4.2 << " : " <<' ' << eval_interval<inter1_t>(4.2) << ' ' << eval_interval<inter2_t>(4.2) << ' ' << eval_interval<merged_t>(4.2) << std::endl;
	std::cout << 9.8 << " : " <<' ' << eval_interval<inter1_t>(9.8) << ' ' << eval_interval<inter2_t>(9.8) << ' ' << eval_interval<merged_t>(9.8) << std::endl;
	std::cout << std::endl;

	return 0;
}

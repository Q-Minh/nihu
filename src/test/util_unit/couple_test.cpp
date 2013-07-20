#include "util/newcouple.hpp"
#include <iostream>
#include <Eigen/Dense>

int main(void)
{
	{
		couple<int, char> a(1, 'b');
		std::cout << a << std::endl;
	}
	
	{
		auto c = create_couple(1, 2);	// by value
		std::cout << c << '\n';
		c.get_first() += 1;
		std::cout << c << '\n';
	}

	{
		int a = 1, b = 2;
		auto c = create_couple(a, b);	// by reference
		std::cout << c << '\n';
		c.get_first() += 1;
		std::cout << c << '\n';
		std::cout << a << '\n';
	}

	{
		int const a = 1;
		int const b = 2;
		auto c = create_couple(a, b);	// by const reference
		std::cout << c << '\n';
	}


	{
		int const a = 1;
		auto c = create_couple(a, 1);	// by const reference
		std::cout << c << '\n';
	}

	{
		auto  mat = Eigen::Matrix<int, 3, 1>::Ones();
		std::cout << (create_couple(1, 2) * mat).get_first() << std::endl;
	}


	{
		auto const c = create_couple(1, 2);
		auto  mat = Eigen::Matrix<int, 3, 1>::Ones();
		std::cout << mat * c * (3 * mat.transpose()) << std::endl;
	}
}


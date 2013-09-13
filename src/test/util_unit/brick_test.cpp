// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "util/brick.hpp"
#include <iostream>
#include <type_traits>

template <int N>
struct n
{
	template <class wall>
	class brick : public wall
	{
	public:
		brick(void) : wall()
		{
			std::cout << "<" << N;
		}

		int get_result(void) const
		{
			return N;
		}
	};
};

typedef glue<n<3>::brick, glue<n<2>::brick, glue<n<1>::brick> > > t1;
typedef glue<n<6>::brick, glue<n<5>::brick, glue<n<4>::brick> > > t2;
typedef glue<n<7>::brick, glue<n<5>::brick, glue<n<0>::brick> > > t3;
typedef merge<t1, t2>::type t1t2;
typedef merge<t1t2, t3>::type t1t2t3;

typedef build<n<1>, n<2>, n<3> >::type t11;
typedef build<n<4>, n<5>, n<6> >::type t22;
typedef build<n<0>, n<5>, n<7> >::type t33;
typedef merge<t11, t22>::type t11t22;
typedef merge<t11t22, t33>::type t11t22t33;

int main(void)
{
	t1 _t1; std::cout << std::endl;
	t2 _t2; std::cout << std::endl;
	t3 _t3; std::cout << std::endl;
	t1t2 _t1t2; std::cout << std::endl;
	t1t2t3 _t1t2t3; std::cout << std::endl;

	std::cout << static_cast<find_in_wall<t1, t1t2t3>::type const &>(_t1t2t3).get_result() << std::endl;
	std::cout << static_cast<find_in_wall<t2, t1t2t3>::type const &>(_t1t2t3).get_result() << std::endl;
	std::cout << static_cast<find_in_wall<build<n<2> >::type, t1t2t3>::type const &>(_t1t2t3).get_result() << std::endl;

	return 0;
}


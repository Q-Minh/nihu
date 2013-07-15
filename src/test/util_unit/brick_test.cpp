#include "util/brick.hpp"
#include <iostream>

template <int N>
struct n
{
	template <class wall>
	class brick : public wall
	{
	public:
		brick(void) : wall()
		{
			std::cout << N << ' ';
		}
	};
};

typedef build<n<3>::brick, build<n<2>::brick, build<n<1>::brick	> > > t1;
typedef build<n<6>::brick, build<n<5>::brick, build<n<4>::brick	> > > t2;
typedef build<n<7>::brick, build<n<5>::brick, build<n<0>::brick	> > > t3;
typedef merge<t1, t2>::type t1t2;
typedef merge<t1t2, t3>::type t1t2t3;


int main(void)
{
	t1 _t1; std::cout << std::endl;
	t2 _t2; std::cout << std::endl;
	t3 _t3; std::cout << std::endl;
	t1t2 _t1t2; std::cout << std::endl;
	t1t2t3 _t1t2t3; std::cout << std::endl;
	
	return 0;
}


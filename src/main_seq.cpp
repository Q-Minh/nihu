#include <iostream>

#include "sequence.hpp"

int main(void)
{
	typedef arg_seq4<int_<1>, int_<2>, int_<100>, int_<-5> > tomb;
	typedef begin<tomb>::type beg;

	return deref<next<next<next<beg>::type>::type>::type>::type::value;
}


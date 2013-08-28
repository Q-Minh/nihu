#include "util/collection.hpp"
#include <iostream>

class empty
{
};

class A
{
public:
	A(int data) : m_data(data) { }
	int get_data(void) const { return m_data; }
private:
	int m_data;
};

class B
{
public:
	B(double data) : m_data(data) { }
	double get_data(void) const { return m_data; }
private:
	double m_data;
};

int main(void)
{
	typedef collect<empty> empty_collection_t;
	typedef collect<A, empty> A_collection_t;
	typedef collect<B> B_collection_t;

	empty_collection_t ec;
	A_collection_t ac(A(3), empty());
	B_collection_t bc(B(2.0));

	std::cout << ac.A::get_data() << std::endl;

	auto m = merge_data(ac, bc);
	std::cout << m.A::get_data() << std::endl;
}

#include <iostream>
using namespace std;

template <class Derived>
class crtp_base
{
public:
	Derived const &derived() const
	{
		return static_cast<Derived const &>(*this);
	}
};

template <class T>
class Aimpl
{
public:
	void implementation(void) const
	{
		cout << "a\n";
	}
};

template <class T>
class A :
	public Aimpl<T>,
	public crtp_base<A<T> >
{
};

class B :
	public crtp_base<B>,
	public Aimpl<char>
{
public:
	void implementation(void) const
	{
		cout << "b\n";
	}
};

template <class T>
void fun(crtp_base<T> const &arg)
{
	arg.derived().implementation();
}


int main(void)
{
	A<char> a;
	B const &bref = static_cast<B const &>(static_cast<Aimpl<char> const &>(a));

	a.implementation();
	bref.implementation();

	fun(a);
	fun(bref);

	return 0;
}


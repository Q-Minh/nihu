#include <type_traits>
#include <iostream>
using namespace std;

// simple model of couple
template <class T>
class storer
{
public:
	storer (T &&t) : m_t(forward<T>(t))
	{
		if (is_lvalue_reference<T>::value)
			cout << "stored by lvalue ref\n";
		else if (is_rvalue_reference<T>::value)
			cout << "stored by rvalue ref\n";
		else
			cout << "stored by value\n";
	}
	
private:
	T m_t;
};

// Gecijó factory. El sem hiszem.
template <class T>
storer<T> factory(T &&t)
{
	return storer<T>(forward<T>(t));
}


// test

struct C
{
	C() : state(true) { cout << "ctor\n"; }
	C(C const &other) : state(other.state) { cout << "copy\n"; }
	C(C &&other) : state(other.state) { other.state = false; cout << "move\n"; }
	~C() { cout << "dtor: " << state << '\n'; }
	bool state;
};

C func(void)
{
	return C();
}

int main(void)
{
	C a;				// instance (ctor)
	C const ca = a;		// instance (copy ctor)
	C &ra = a;			// reference
	C const &cra = a;	// const reference
	C && rra = func();	// rvalue reference to temporary (ctor)

	cout << "\n\n";
	
	factory(ca); cout << "\n";
	factory(cra); cout << "\n";
	
	factory(C()); cout << "\n";
	factory(func()); cout << "\n";
	factory(forward<C>(rra)); cout << "\n";
	
	factory(a); cout << "\n";
	factory(ra); cout << "\n";

	return 0;
}


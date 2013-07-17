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

// converts to T from T&&, otherwise unchanged. This is the store type of the couple
template <class T>
struct remove_rvalue_reference
{
	typedef T type;
};

template <class T>
struct remove_rvalue_reference<T&&>
{
	typedef T type;
};


// Gecijó factory. El sem hiszem.
template <class T>
storer<typename remove_rvalue_reference<T&&>::type>
	factory(T &&t)
{
	return storer<typename remove_rvalue_reference<T&&>::type>(forward<T>(t));
}


// test

struct C
{
	C(C const &other) { cout << "copy\n"; }
	C(C &&other) { cout << "move\n"; }
	C() { cout << "ctor\n"; }
	~C() { cout << "dtor\n"; }
};

C func(void)
{
	return C();
}

int main(void)
{
	C a;
	C const ca = a;
	C &ra = a;
	C const &cra = a;
	
	cout << "\n\n";
	
	factory(ca); cout << "\n";
	factory(cra); cout << "\n";
	
	factory(C()); cout << "\n";
	factory(func()); cout << "\n";
	
	factory(a); cout << "\n";
	factory(ra); cout << "\n";

	return 0;
}


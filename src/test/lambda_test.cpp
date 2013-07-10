#include "../tmp/lambda.hpp"
#include "../tmp/integer.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"
#include <iostream>
#include <typeinfo>

using namespace tmp;
using namespace std;

template <class A, class B>
struct Fuso
{
	typedef A type;
};

template <class A>
struct Transform1 { struct type
{
	void operator() ()
	{
		cout << "doing for input " << A::value << endl;
	}
};};


template <class A, class B>
struct Transform { struct type
{
	void operator() ()
	{
		cout << "doing for inputs " << A::value << " " << B::value << endl;
	}
};};

int main(void)
{
/*
	// test isPlaceholder
	static_assert(isPlaceholder<_1>::type::value == true, "isPlaceholder error");
	static_assert(isPlaceholder<int>::type::value == false, "isPlaceholder error");

	// test isPlaceholderExpression
	static_assert(isPlaceholderExpression<_1>::type::value == true, "isPlaceholderExpression error");
	static_assert(isPlaceholderExpression<int>::type::value == false, "isPlaceholderExpression error");
	static_assert(isPlaceholderExpression<Fuso<int, char> >::type::value == false, "isPlaceholderExpression error");
	static_assert(isPlaceholderExpression<Fuso<_1, char> >::type::value == true, "isPlaceholderExpression error");
	static_assert(isPlaceholderExpression<Fuso<int, _1> >::type::value == true, "isPlaceholderExpression error");
	static_assert(isPlaceholderExpression<Fuso<_2, _1> >::type::value == true, "isPlaceholderExpression error");

	// test lambda
	typedef lambda<Fuso<int_<20>, _1> >::type metaFunctr1;
	static_assert(metaFunctr1::apply<int_<1> >::type::value == 20, "lambda error");
	typedef lambda<Fuso<_1, int_<20> > >::type metaFunctr2;
	static_assert(metaFunctr2::apply<int_<1> >::type::value == 1, "lambda error");
	typedef lambda<Fuso<_1, _2> >::type metaFunctr3;
	static_assert(metaFunctr3::apply<int_<28>, int_<2> >::type::value == 28, "lambda error");
	static_assert(is_same<lambda<int>::type, int>::type::value, "lambda error");

	// test apply
	static_assert(apply<Fuso<int_<20>, _1>, int_<1> >::type::value == 20, "apply error");
	static_assert(apply<Fuso<_1, int_<20> >, int_<1> >::type::value == 1, "apply error");
	static_assert(apply<Fuso<_1, _2>, int_<1>, int_<2> >::type::value == 1, "apply error");
	

	typedef vector<int_<1>, int_<2>, int_<3> > v1;
	typedef vector<int_<11>, int_<12>, int_<13> > v2;
*/


/*
	// test call_each
	call_each<v1, Transform<_1, int_<0> > >();
*/

	cout << endl;


/*
	typedef lambda<Transform1<_1 > >::type meta1functor;
	typedef meta1functor::apply<int_<1> > partial1;
	static_assert(isPlaceholderExpression<partial1>::type::value == false, "partial1 is plhexpr");

	partial1 p1;
	cout << typeid(p1).name() << endl;
*/

	typedef lambda<Transform<_1, _2> >::type metafunctor;

	typedef internal::lambda_plExp<metafunctor::apply<int_<1>, _1> >::type l;

/*
	typedef l::apply<int_<2> >::type cur;
	static_assert(isPlaceholderExpression<cur>::type::value == false, "cur is plhexpr");
	cur c;
	cout << typeid(c).name() << endl;
*/

/*

	call_each<v2, partially_evaluated_transform, char, int>('a', 1);
*/

	


/*
	// test d_call_each
	d_call_each<v1, v2, Transform<_1,_2>, char, int>('a', 1);
*/
	return 0;
}


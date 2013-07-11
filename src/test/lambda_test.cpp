#include "../tmp/lambda.hpp"
#include "../tmp/integer.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"
#include <iostream>
#include <typeinfo>

using namespace tmp;
using namespace std;

template <class A>
struct MF1 {
	typedef A type;
};

template <class A, class B>
struct MF2 {
	typedef B type;
};

template <class A>
struct TR1 
{ 
	struct type
	{
		void operator() ()
		{
			cout << "doing for input " << A::value << endl;
		}
	};
};


template <class A, class B>
struct TR2 
{ 
	struct type
	{
		void operator() ()
		{
			cout << "doing for inputs " << A::value << " " << B::value << endl;
		}
	};
};

int main(void)
{
	// Predeclaration of vector types
	typedef vector<int_<1>, int_<2>, int_<3> > VEC1;
	typedef vector<int_<-11>, int_<-23>, int_<-38> > VEC2;
	
	// One argument test
	typedef apply<TR1<_1 >, int_<2> >::type MF1_APP;
	MF1_APP mf1_app;
	mf1_app();
	
	
	call_each<VEC1, TR1<_1 > >();
	
	// Two arguments test
	typedef apply<TR2<_1, int_<-32> >, int_<2> >::type MF2_APP;
	MF2_APP mf2_app;
	mf2_app();
	
	
	call_each<VEC1, TR2<_1, int_<-37> > >();
	typedef apply<TR2 <_1, _2>, int_<-69>, int_<-879> > ::type MF2_FULL;
	MF2_FULL mf2_full;
	mf2_full();
	
	typedef internal::lambda_plExp<TR2<_1, _2>>::type test;
	typedef test::apply<int_<-88>, int_<-99>>::type apply_test;
	apply_test at;
	at();
	
	
	typedef vector<int_<3> > VEC3;
	typedef vector<int_<4> > VEC4;

	d_call_each<VEC1, VEC2, TR2<_1, _2 > >();
	
	//typedef apply<test::apply<int_<77> >, int_<-3> > test2;
	//test2 t2;
	//t2();
	

	
/*
	typedef lambda<Transform1<_1 > >::type meta1functor;
	typedef meta1functor::apply<int_<1> > partial1;
	static_assert(isPlaceholderExpression<partial1>::type::value == false, "partial1 is plhexpr");

	partial1 p1;
	cout << typeid(p1).name() << endl;
*/
/*
	typedef lambda<Transform<_1, _2> >::type metafunctor;

	typedef internal::lambda_plExp<metafunctor::apply<int_<1>, _1> >::type l;
*/
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


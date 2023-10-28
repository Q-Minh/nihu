/**
 * \file metaprogramming.hpp
 * \brief Introduces the metaprogramming concepts of NiHu
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

//! [MetaFunction]
struct AMetaFunction
{
	typedef int type;
};
//! [MetaFunction]

//! [MetaFunctionClass]
struct AMetaFunctionClass
{
	template <typename T>
	struct apply
	{
		typedef T type;
	};
};
//! [MetaFunctionClass]

//! [ArgumentSelector]
template <unsigned N, class T, class ...Args>
struct select_argument : select_argument<N-1, Args...> {};

template <class T, class ...Args>
struct select_argument<1U, T, Args...>
{
	typedef T type;
};
//! [ArgumentSelector]

//! [PlaceHolder]
template <unsigned N>
struct arg
{
	template <class...Args>
	struct apply : select_argument<N, Args...> {};
};
//! [PlaceHolder]
	
//! [PlaceHolderTypedef]
typedef arg<1> _1;
typedef arg<2> _2;
...
//! [PlaceHolderTypedef]

//! [PlaceHolderDecide]
template <class C>
struct isPlaceholderExpression : std::false_type {};

template <unsigned N>
struct isPlaceholderExpression<arg<N> > : std::true_type {};

template <template <class...Args> class MF, class...Args>
struct isPlaceholderExpression< MF<Args...> > :
	   containsPlaceholderExpression<Args...> {};
//! [PlaceHolderDecide]

//! [LambdaExample]	   
template <class T>	   
struct MyMF {
	struct type {
		// Assuming T defines value as an integer expression
		void operator() (void) {cout << T::value << endl;}
	};
};

typedef vector<int_<1>, int_<2>, int_<3> > vec; 	// A type vector
typedef MyMF<_1> lazy_exp;							// Type expression with lazy evaluation
call_each<vec, lazy_exp >();						// Evaluation performed here
//! [LambdaExample]
	   
//! [LambdaMetaFun]
template <class Exp>
	struct lambda : if_<
		typename isPlaceholderExpression<Exp>::type,
		typename internal::lambda_plExp<Exp>::type,
		Exp
	> {};
//! [LambdaMetaFun]
	
//! [LambdaMetaFunExample]
typedef MyMF<_1> lazy_exp;							// Type expression with lazy evaluation
typedef lambda<lazy_exp>::type lambda_t;			// Wrap as lambda metafunction
typedef lambda_t::apply<int_<4> >::type eval_t;		// The evaluated type
eval_t e; 											// Instantiate 
e();												// Call evaluation operator
//! [LambdaMetaFunExample]
	   
//! [ApplyMetaFunction]
template <class Fun, class...Args>
struct apply : lambda<Fun>::type::template apply<Args...> {};
//! [ApplyMetaFunction]

//! [TypedefForwarding]
struct parent_struct {
	typedef int type;
}

struct child_struct : parent_struct {};
//! [TypedefForwarding]

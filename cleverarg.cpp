#include <iostream>
#include <type_traits>

#include "src/tmp/bool.hpp"

// a simple kernelinput
class ArgA
{
public:
	// constructor from int
	ArgA(int par) :
		m_a(par)
	{
		std::cout << "ArgA ctor\n";
	}
	
	// return stored value
	int get_a(void) const
	{
		return m_a;
	}

private:
	int m_a;
};

// A kernel depending on the simple input
class KernelA
{
public:
	// the input type
	typedef ArgA arg_t;
	
	// evaluation
	void eval(arg_t const &arg)
	{
		m_value = arg.get_a();
		std::cout << "KernelA evaluated\n";
	}
	
	// return stored kernel value
	int get_value(void) const
	{
		return m_value;
	}

private:
	int m_value;
};

// a complex kernelinput derived from the simple one
class ArgB :
	public ArgA
{
public:
	// constructor from int
	ArgB(int par) :
		ArgA(par), m_b(2*par)
	{
		std::cout << "ArgB ctor\n";
	}
	
	// return stored value
	int get_b(void) const
	{
		return m_b;
	}

private:
	int m_b;
};

// A complex kernel depending on the complex input and on the simple kernel
class KernelB :
	public KernelA
{
public:
	// the input type
	typedef ArgB arg_t;

	// evaluate using KernelA
	void eval(arg_t const &arg)
	{
		KernelA::eval(arg);
		m_value = arg.get_b() * KernelA::get_value();
		std::cout << "KernelB evaluated\n";
	}
	
	// return stored result
	int get_value(void) const
	{
		return m_value;
	}

private:
	int m_value;
};


// metafunction returning a combined input from A and B with a constructor
template <class A, class B, class...CtorArgs>
struct combined_input
{
	// return a class derived from both A and B
	struct type :
		public A, public B
	{
		// constructor
		type(CtorArgs...args) :
			A(args...), B(args...)
		{
		}
	};
};

// specialisation of combined_input for the trivial case A = B
template <class A, class...CtorArgs>
struct combined_input<A, A, CtorArgs...>
{
	// return A
	typedef A type;
};


// return smallest common wrapper input for inputs A and B
template <class A, class B, class CtorArg>
struct smallest_input_wrapper : tmp::if_<
	typename std::is_base_of<A, B>::type, B, typename tmp::if_<
		typename std::is_base_of<B, A>::type, A, typename combined_input<A, B, CtorArg>::type
	>::type
> {};



// metafunction returning a combined kernel from A and B with a constructor and a typedef
template <class A, class B, class InputArg>
struct combined_kernel
{
	// return a class derived from both A and B
	struct type :
		public A, public B
	{
		// the eval argument type
		typedef InputArg arg_t;
		
		// constructor
		type(A const &a, B const &b) :
			A(a), B(b)
		{
		}
		
		// independent sequential evaluations
		void eval(InputArg const &arg)
		{
			A::eval(arg);
			B::eval(arg);
		}
	};
};


// specialisation of combined_kernel for the A = B case
template <class A, class InputArg>
struct combined_kernel<A, A, InputArg>
{
	// return a class derived from only A
	struct type :
		public A
	{
		// the eval argument type
		typedef InputArg arg_t;
		
		// constructor
		type(A const &a, A const &) :
			A(a)
		{
		}

		using A::eval;
	};
};


// specialisation of combined_kernel for the A includes B case
template <class A, class B, class InputArg>
struct Asuperior
{
	// return a class derived from only A
	struct type :
		public A
	{
		// the eval argument type
		typedef InputArg arg_t;
		
		// constructor
		type(A const &a, B const &) :
			A(a)
		{
		}
		
		using A::eval;
	};
};


// specialisation of combined_kernel for the B includes A case
template <class A, class B, class InputArg>
struct Bsuperior
{
	// return a class derived from only B
	struct type :
		public B
	{
		// the eval argument type
		typedef InputArg arg_t;
		
		// constructor
		type(A const &, B const &b) :
			B(b)
		{
		}
		
		using B::eval;
	};
};


// return smallest common wrapper kernel for kernels A and B
template <class A, class B, class InputArg>
struct smallest_kernel_wrapper : tmp::if_<
	typename std::is_base_of<A, B>::type,
	typename Bsuperior<A, B, InputArg>::type,
	typename tmp::if_<
		typename std::is_base_of<B, A>::type,
		typename Asuperior<A, B, InputArg>::type,
		typename combined_kernel<A, B, InputArg>::type
	>::type
> {};


// combine two kernels, provide an optimal kernel input
template <class K1, class K2>
struct super_kernel
{
	// the common kernel input
	typedef typename smallest_input_wrapper<
		typename K1::arg_t,
		typename K2::arg_t,
		int
	>::type argt;

	// the common kernel
	typedef typename smallest_kernel_wrapper<
		K1,
		K2,
		argt const &	// we use the common kernel input here
	>::type type;
};


// superkernel factory operator
template <class KA, class KB>
typename super_kernel<KA, KB>::type
	operator,(KA const &ka, KB const &kb)
{
	return typename super_kernel<KA, KB>::type(ka, kb);
}


template<class K>
void test(K k)
{
	typename K::arg_t arg(2);
	k.eval(arg);
}

int main(void)
{
	KernelA ka;
	KernelB kb;
	
	std::cout << std::endl << "A-B\n";
	test( (ka, kb) );
	
	std::cout << std::endl << "B-A\n";
	test( (kb, ka) );
	
	std::cout << std::endl << "A-A\n";
	test( (ka, ka) );
	
	std::cout << std::endl << "B-B\n";
	test( (kb, kb) );
	
	return 0;
}


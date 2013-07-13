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

		// eval only evaluates A
		void eval(InputArg const &arg)
		{
			A::eval(arg);
		}
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


// return smallest common wrapper input for inputs A and B
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


// return smallest common wrapper input for inputs A and B
template <class A, class InputArg>
struct smallest_kernel_wrapper<A, A, InputArg> :
	combined_kernel<A, A, InputArg> {};


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


// superkernel factory function
template <class KA, class KB>
typename super_kernel<KA, KB>::type
	create_super_kernel(KA const &ka, KB const &kb)
{
	return typename super_kernel<KA, KB>::type(ka, kb);
}



int main(void)
{
	// instantiate two separate kernels (in main module with wave numbers)
	KernelA ka;
	KernelB kb;
	
	// create the super kernel instance (still in main module)
	auto kc = create_super_kernel(kb, ka);
	
	// retreive optimal input type from superkernel and instantiate input (in integration module)
	decltype(kc)::arg_t arg(2);

	// evaluate superkernel with the optimal argument
	kc.eval(arg);

	// get return value of first kernel
	std::cout << kc.KernelA::get_value() << std::endl;
	// get return value of second kernel
	std::cout << kc.KernelB::get_value() << std::endl;
	
	// GEEEEXCIIIIII !!!

	return 0;
}


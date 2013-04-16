Homogeneous collections {#homocollect}
=======================

[TOC]

Polymorphism in NiHu {#nihupoly}
====================

Evaluating the Weighted Residual approach means integrating a kernel over a set of elements extended with some shape functions. Obviously, when programming a BEM, we want to write our integration routine generally, so that our code remains capable to handle as many different element and kernel types as possible.

Dynamic polymorphism {#dynpoly}
--------------------

A straightforward C++ solution to this problem is *dynamic polymorphism* implemented with virtual functions, abstract base classes and heterogeneous collections. Using this technique, an integration routine can be written that receives pointers to abstract kernel and element interface classes as input, and performs integration by invoking the functionalities of the abstract interfaces. The specific behaviour of different subproblems are implemented in specialised derived classes of the abstract base.

The main advantage of dynamic polymorphism is that it can easily and transparently manage heterogeneous problmes, for example BEM problems where the mesh contains different kind of elements. However, dynamic polymorphism has some important drawbacks too:
- The actual implementation functions are dispatched during runtime, so function calls cannot be inlined by the optimising compiler.
- Although different specialisations of the same abstract base can differ in implementation, they cannot differ in argument types and return types. This means that if we would like to keep our implementation general, the interface data needs to be expressed in a general way, which usually results in hard-to-optimize dynamically allocated or redundant structures.

These drawbacks can add up to a significant decrease of overall performance.

Static polymorphism {#statpoly}
--------------------

_Static polymorphism_, on the contrary, means that we implement our problem generally, but keep in mind that only one single specific instance of the general problem will be compiled and run at once. If we want to invoke different special instances of the same abstract problem, we need to compile and call the different specialisations separately.

In C++ static polymorphism is generally implemented using class templates and function templates. In the example below we demonstrate its application with a simple class template that integrates an arbitrary kernel over an arbitrary element:

~~~~~~~~{.cpp}
template <class Elem, class Kernel>
class weighted_integral
{
public:
	// compute the integral's result type
	typedef typename product_type<Kernel::result_t, Elem::shape_function_t>::type result_t;

	// integrate a kernel over an element
	static result_t integrate(Elem const &elem, Kernel const &kernel)
	{
		result_t result = result_t();	// clear the result
		// traverse the element's quadrature points
		for (auto it = elem.quadr_begin(); it != elem.quadr_end(); ++it)
		{
			result +=
				kernel.eval(*it) *	// evaluate kernel
				elem.get_shape(*it);	// evaluate element shape functions
		}
	}
	return result;
};
~~~~~~~~

When reading the above code, the following observations should be made:
- Class `weighted_resudial` is implemented in terms of a general `Kernel` and `Elem` class, but it is going to be instantiated for a single specific `Kernel` and `Elem` type.
- The `integrate` function's arguments and return type are specific to the template parameters `Elem` and `Kernel`. This means that different specialisations can have different interface, and the function can be inlined and optimised.

The class template can be instantiated for any `Elem` and `Kernel` type (that provide an appropriate `eval` function, `get_shape` function and quadrature iterators). For example, specialisation with a Helmholtz kernel and triangle elements can be invoked by compiling the code segment below:

~~~~~~~~{.cpp}
helmholtz_kernel kernel;

std::vector<tria_elem> tria_elements;

for (auto it = tria_elements.begin(); it != tria_elements.end(); ++it)
	weighted_integral<tria_elem, helmholtz_kernel>::integrate(*it, kernel);
~~~~~~~~

What if we need to incorporate quad elements into our program? We need to rewrite the code as

~~~~~~~~{.cpp}
helmholtz_kernel kernel;

std::vector<tria_elem> tria_elements;
std::vector<quad_elem> quad_elements;

for (auto it = tria_elements.begin(); it != tria_elements.end(); ++it)
	weighted_integral<tria_elem, helmholtz_kernel>::integrate(*it, kernel);
for (auto it = quad_elements.begin(); it != quad_elements.end(); ++it)
	weighted_integral<quad_elem, helmholtz_kernel>::integrate(*it, kernel);
~~~~~~~~
It is obvious that integrating tria and quad elements can only be done separately, one element type after the other. This is the price paid for optimal performance.

What if we need to evaluate double integrals (Galerkin BEM) with arbitrary test and trial shape functions, we need to handle three different kernels and five different element types? We have to repeat our simple traversing code segment \f$5^2\times3=75\f$ times! This is the point where template metaprogramming needs to be exploited.

Template Metaprogramming {#tmp}
------------------------

Metaprogramming is writing programs that write programs. Template Metaprogramming is using C++ templates to write programs that write programs. To generalise the above example, we need two things:
1. A unified container that contains a vector of tria elements, a vector of quad elements and further vectors of other element types.
2. A loop over types that can be called to traverse each container, instantiate class template `weighted_integral` for each element type and invoke the `integrate` function for each element of the subcontainer.

### Inheriters {#inheriter}

The first task is accomplished by an inheritance trick:

~~~~~~{.cpp}
class A1
{
protected:
	std::vector<tria_elem> container;
};

class A2 : public A1
{
protected:
	std::vector<quad_elem> container;
};

class A3 : public A2
{
protected:
	std::vector<any_other_elem> container;
};

...
~~~~~~

As each class `Ai` is derived from the previous class `A(i-1)`, class `A3` contains each element type vectors. The general recursive implementation of this inheritance pattern can be written as

~~~~~~{.cpp}
template <class TypeVector>
class Container : public Container<typename tail<TypeVector>::type>
{
protected:
	std::vector<typename head<TypeVector>::type> container;
};

// terminating condition
template <>
class Container<EmptyType> {};
~~~~~~

This is the technique how NiHu stores inhomogeneous meshes. For more information check out the class documentattion of class ::Mesh and its inhomogeneous container member ::Mesh::m_elements.

### call_each {#calleach}


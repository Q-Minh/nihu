Polymorphism in NiHu {#polymorphism}
====================

[TOC]

[metafunction]:http://www.boost.org/doc/libs/1_53_0/libs/mpl/doc/refmanual/metafunction.html

This page explains a bit about how NiHu works.
The topic discussed here is not about Boundary Elements but programming.
We explain how static polymorphism and c++ template metaprogramming is applied to efficiently incorporate inhomogeneous meshes into the NiHu toolbox.

Polymorphism in NiHu {#nihupoly}
====================

Evaluating the Weighted Residual approach means integrating a kernel \f$K(x_0,x)\f$ over a set of elements \f$E\f$ extended with some shape functions \f$N(x)\f$:
\f[
R = \int_{E}K(x_0,x) N(x) dx
\f]
Obviously, when programming a BEM, we want to write our integration routine generally, so that our code remains capable to handle as many different element and kernel types as possible.

Dynamic polymorphism {#dynpoly}
--------------------

A straightforward C++ solution to this problem is *dynamic polymorphism* implemented with virtual functions, abstract base classes and heterogeneous collections.
Using this technique, the integration routine receives pointers to abstract kernel and element interface classes as input, and performs integration by invoking the functionalities of the abstract interfaces.
The specific behaviour of different subproblems are implemented in specialised derived classes of the abstract base.

The main advantage of dynamic polymorphism is that it can easily and transparently manage heterogeneous problems, for example BEM problems where the mesh contains different kind of elements.
However, dynamic polymorphism has some important drawbacks too:
- The actual implementation functions are dispatched during runtime, so function calls cannot be inlined by the optimising compiler.
- Although different specialisations of the same abstract base can differ in implementation, they cannot differ in argument types and return types. So in order to keep your implementation general, the interface data needs to be expressed in a general way, which usually results in hard-to-optimize dynamically allocated or redundant structures.

These drawbacks can add up to a significant decrease of overall performance.

Static polymorphism {#statpoly}
--------------------

_Static polymorphism_, on the contrary, means that we still implement our problem generally, but keep in mind that only one single specific instance of the general problem will be compiled and run at once.
As dispatch is done compile time, invoking different special instances of the same abstract problem is managed by compiling and calling the different specialisations separately.

In C++ static polymorphism can be implemented using class templates and function templates.
In the example below we demonstrate its application with a simple class template that integrates an arbitrary kernel over an arbitrary element:

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
		return result;
	}
};
~~~~~~~~

When reading the above code, the following observations should be made:
- Class template `weighted_integral` is implemented in terms of a general `Kernel` and `Elem` class, but it is going to be instantiated for a single specific `Kernel` and `Elem` type.
- The `integrate` function's arguments and return type are specific to the template parameters `Elem` and `Kernel`. This means that different specialisations can have different interface, and the function can be inlined and optimised.
- The return type of function `integrate` is "computed" by a [metafunction] `product_type`.

The class template can be instantiated for any `Elem` and `Kernel` type (that provide an appropriate `eval` function, `get_shape` function and quadrature iterators).
For example, specialisation with a Helmholtz kernel and triangle elements can be invoked by compiling the code segment below:

~~~~~~~~{.cpp}
helmholtz_kernel kernel;

std::vector<tria_1_elem> tria_elements;

for (auto it = tria_elements.begin(); it != tria_elements.end(); ++it)
	weighted_integral<tria_1_elem, helmholtz_kernel>::integrate(*it, kernel);
~~~~~~~~

What if we need to incorporate quad elements into our program? We need to rewrite the code as

~~~~~~~~{.cpp}
helmholtz_kernel kernel;

std::vector<tria_1_elem> tria_elements;
std::vector<quad_1_elem> quad_elements;

for (auto it = tria_elements.begin(); it != tria_elements.end(); ++it)
	weighted_integral<tria_1_elem, helmholtz_kernel>::integrate(*it, kernel);
for (auto it = quad_elements.begin(); it != quad_elements.end(); ++it)
	weighted_integral<quad_1_elem, helmholtz_kernel>::integrate(*it, kernel);
~~~~~~~~
It is obvious that integrating tria and quad elements can only be done separately, one element type after the other. This is the price paid for better performance.

What if we need to evaluate double integrals (Galerkin BEM) with arbitrary test and trial shape functions, five different element types and three different kernels?
We have to repeat our simple traversing code segment \f$5^2\times3=75\f$ times!
This is the point where template metaprogramming becomes inevitable.

Template Metaprogramming {#tmp}
------------------------

Metaprogramming is writing programs that write programs.
Template Metaprogramming is using C++ templates to write programs that write programs.
In the following we will apply template metaprogramming to code generation.

To generalise the above example, we need two things:
1. A unified container that contains a `std::vector<tria_1_elem>` member, a `std::vector<quad_1_elem>` member and further vectors of other element types.
2. A loop over types that can be invoked to traverse each container, instantiate class template `weighted_integral` for each element type and call the `integrate` function for each element of the element type's vector container.

### Inheriters {#inheriter}

The first task is accomplished by an inheritance trick.
Our container class has to be derived from all the different `std::vector<>` containers:

~~~~~~{.cpp}
class A : public std::vector<tria_1_elem>, std::vector<quad_1_elem>, std::vector<any_other_elem> {};
~~~~~~

The above pattern is easily implemented using a series of classes inheriting from at most two bases
~~~~~~{.cpp}
class C : public std::vector<any_other_elem> {};
class B : public std::vector<quad_1_elem>, C {};
class A : public std::vector<tria_1_elem>, B {};
~~~~~~

The recursive implementation of this pattern is the following:
~~~~~~{.cpp}
template <class TypeVector>
class Container :
	public std::vector<typename head<TypeVector>::type>,
	Container<typename tail<TypeVector>::type>
{};

// terminating condition
template <>
class Container<EmptyType> {};
~~~~~~

And the instantiation for a given vector of element types is as follows:
~~~~~~{.cpp}
// the element type vector
typedef tmp::vector<tria_1_elem, quad_1_elem, any_other_elem> elem_vector_t;
// the inheriter container
typedef Container<elem_vector_t> mesh_t;
// the mesh instance
mesh_t mesh;
~~~~~~

This is the technique how NiHu stores inhomogeneous meshes.
For more information check out the class documentattion of class ::Mesh and its inhomogeneous container member ::Mesh::m_elements.

### call_each {#calleach}

Now that we have our inhomogeneous container `mesh` consisting of homogeneous `std::vector`-s, we can rewrite our code segment that integrates over all the element types

~~~~~~{.cpp}
for (auto it = mesh.std::vector<tria_1_elem>::begin(); it != mesh.std::vector<tria_1_elem>::end(); ++it)
	weighted_integral<tria_1_elem, helmholtz_kernel>::integrate(*it, kernel);
for (auto it = mesh.std::vector<quad_1_elem>::begin(); it != mesh.std::vector<quad_1_elem>::end(); ++it)
	weighted_integral<quad_1_elem, helmholtz_kernel>::integrate(*it, kernel);
for (auto it = mesh.std::vector<any_other_elem>::begin(); it != mesh.std::vector<any_other_elem>::end(); ++it)
	weighted_integral<any_other_elem, helmholtz_kernel>::integrate(*it, kernel);
~~~~~~

This pattern can be generalised using NiHu's `tmp::call_each` code generating metafunction.
`call_each` is used to instantiate a class template for each entry of a specified type vector, and call the instantiated class' nested functor called `type`.
The application is demonstrated as follows:

~~~~~~{.cpp}
// the class template
template <class elem_t>
struct integrator
{
	// the nested functor
	struct type
	{
		void operator() (kernel_t kernel, mesh_t &mesh)
		{
			for (auto it = mesh.std::vector<elem_t>::begin(); it != mesh.std::vector<elem_t>::end(); ++it)
				weighted_integral<elem_t, kernel_t>::integrate(*it, kernel);
		}
	};
};

tmp::call_each<
	tmp::vector<tria_1_elem, quad_1_elem, any_other_elem>,	// the element type vector
	integrator<_1>,	// the class template with a placeholder
	kernel_t,		// the first argument type of the functor
	mesh_t &		// the second argument type of the functor
>(kernel, mesh);	// the arguments
~~~~~~

Homogeneous collections {#homocollect}
=======================

[TOC]

Polymorphism in NiHu {#nihupoly}
====================

Evaluating the Weighted Residual approach means integrating a kernel over a set of elements extended with some shape functions. Obviously, when programming a BEM, we want to write our integration routine generally, so that our code remains capable to handle as many different element and kernel types as possible.

Dynamic polymorphism {#dynpoly}
--------------------

A straightforward C++ solution to this problem is *dynamic polymorphism* implemented with virtual functions. Using this technique, we write an integration routine that receives a pointer to an abstract kernel interface class and a pointer to an abstract element interface class as input argument. The interface classes only declare the functionalities needed to evaluate the integrals, and their derived specialisations implement the actual algorithms. The integration routine invokes the functionalities of the abstract interfaces, and the actual implementations are searched during runtime from virtual function tables.

~~~~~~{.cpp}
result_matrix integrate(Kernel *kernel, Element *elem)
{
	result_matrix res = 0;
	for (it = elem->begin(); it != elem->end(); ++it)
		res += kernel->eval(*it) * elem->shape_function(*it);
	return res;
}
~~~~~~

The main advantage of dynamic polymorphism is that if we build a vector of element pointers, integration over the whole mesh can be performed with a single loop over our container.

~~~~~~{.cpp}
Element *mesh[100];
for (size_t i = 0; i < 100; ++i)
	integrate(kernel, mesh[i]);
~~~~~~

The main drawbacks of the dynamic polymorphism approach are the following:
- The actual function calls are dispatched during runtime, so they cannot be inlined by the compiler.
- The size of the result and temporary vectors is unknown compile time, so they need to be stored dynamically, and even simple traversing loops cannot be unrolled by the compiler.

These drawbacks can add up to a significant decrease of overall performance.


Static polymorphism {#statpoly}
--------------------

Static polymorphism means that we implement our integration routine generally
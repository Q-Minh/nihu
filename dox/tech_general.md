General Techniques {#tech_general}
==================

\page tech_general General Techniques

[TOC]

Store pattern {#tech_general_store}
=============

The store pattern is used to automatically create static instances of class template specialisations.
The pattern is simply implemented as follows:
~~~~~~~~~~~
template <class C>
struct store
{
	static const C m_data;
};

template <class C>
const C store<C>::m_data;
~~~~~~~~~~~

The store pattern's usage is explained with a simple example.
The following code snipped defines a cache template that stores 1000 values of an integral type in a dynamically allocated array.
The cache class is indexable using the overloaded index operator.

\snippet store_test.cpp cache

The following main function uses the store pattern to get two elements from `cache<int>` and one from `cache<char>` This is simply accomplished by instantiating the `store` template with `cache<int>` and `cache<char>` and use their static member `m_data` as follows:

\snippet store_test.cpp main

The code's output is the following:

	cache constructor i
	cache constructor c
	5
	25
	!
	cache destructor c
	cache destructor i

Apparently, both cache's are automatically constructed _once_ at the beginning of the program and are destructed at the end.
This property is useful when different cache's are extensively used throughout different segments of a compilation unit.
The programmer does not need to define the static cache instances, the necessary constructors and destructors are called automatically and only once at the beginning and the end of the program.

This technique is extensively used in NiHu to conveniently define accelerator caches containing quadratures and precomputed weighting functions.
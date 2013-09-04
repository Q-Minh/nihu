C++11 techniques in NiHu {#tech_cpp11}
========================

\page tech_cpp11

[TOC]

Introduction {#tech_cpp11_intro}
============

This document presents the C++11 techniques used in NiHu.

Variadic templates {#tech_cpp11_variadic}
==================

A new feature of C++11 enables the usage of template classes with variable number of template arguments, a.k.a. _variadic templates_.
Variadic templates are created using the ellipses (`...`) syntax as demonstrated in the example below.
The example also shows function parametrisation and calling using a variadic list of arguments.

\snippet cpp11techniques.hpp VariadicTemplateClass

To interpret a code that utilises the variadic syntax you can always think of the ellipses as an extraction of the underlying operation that is the operation applied on the elements of the variadic expression, separated by commas.
For example the expression 

\snippet cpp11techniques.hpp VariadicExample

with `Args = bool, char` is interpreted as

\snippet cpp11techniques.hpp VariadicInterpretation

It is also worth mentoning that the expression `template <class...Args>` also covers the empty `Args` case, which means no template parameters in case of the above example class.

What makes variadic templates extremely useful is, of course, specialisation.
The basic example defines the disjunction of boolean template parameters as

\snippet cpp11techniques.hpp VariadicSpec1
	
As it is seen, the default value is `false`, which will be utilised in the specialisation.
The powerfulness of variadic templates are demonstrated by the specialised classes:

\snippet cpp11techniques.hpp VariadicSpec2

As a result, all the specialisations with any number of arguments such as

\snippet cpp11techniques.hpp VariadicSpec3

will work, and give a correct result.

The auto keyword {#tech_cpp11_auto}
================

C++11 allows to use the `auto` keyword for the definition of variables with explicit initialisation, where the compiler can deduce the type name from the initialisation expression:

\snippet cpp11techniques.hpp Auto first

This new feature becomes really useful when massively templated codes are used. Consider the following example: Class `dirac_wrapper` is used to wrap a function space type into a Dirac function space type. Class `dirac_wrapper` is a template, and `dirac_wrapper` objects can be constructed from any kind of function space objects, as shown in the following snippet.

\snippet cpp11techniques.hpp Dirac wrapper class

Now suppose that we have a variable `func_space` of type `func_space_t`. The corresponding Dirac wrapper object can be constructed as

\snippet cpp11techniques.hpp Dirac wrapper old declaration

This declaration can be made more readable by introducing a factory function template `dirac`:

\snippet cpp11techniques.hpp Dirac wrapper factory

and initialising our new variable with the factory function's return value:

\snippet cpp11techniques.hpp Dirac wrapper usage

note that the compiler deduced the template parameter of the factory function from its argument, and the return type defined the type of the explicit initialiser of the variable definition. The result is a simple Matlab-like readable expression that compiles to the same as the old-fashioned explicit constructor call.

The type_traits library {#tech_cpp11_type_traits}
=======================

lambda expressions {#tech_cpp11_lambda}
==================


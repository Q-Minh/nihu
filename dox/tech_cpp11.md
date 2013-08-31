C++11 techniques in NiHu {#tech_cpp11}
========================

\page tech_cpp11

[TOC]

Introduction {#tech_cpp11_intro}
============

This document presents the C++11 techniques used in NiHu.

Variadic templates {#tech_cpp11_variadic}
==================

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


Metaprogramming concepts in NiHu {#metaprogramming}
================================

[TOC]

Introduction {#intro}
============

This document briefly summarises the metaprogramming concepts used in NiHu. Understanding the code snippets contained in this document requires your being familiar with the usage of variadic templates, introduced in the C++11 standard. If you are not familiar with variadic templates it is advised to read [C++11 techniques in NiHu](cpp11techniques.md) first.

Nomenclature {#nomen}
============

In the following some metaprogramming elements are explained. Most of the nomenclature explained herein is very similar to the definitions of boost, so if you are familiar with boost, you might as well skip this section.

Metafunctions {#metafun}
-------------

Metafunctions are program elements (classes in C++) that can return computed types. Whereas a function can receive its arguments as values and computes a return value, a metafunction can receive its arguments as template parameters and computes a return type. One great difference is that while functions compute their return value in runtime, metafunctions must be evaluated in compile time. Metafunctions are implemented in NiHu usings simple structures with nested type definitions, called `type`.

\snippet metaprogramming.hpp MetaFunction

Argument selector {#argselector}
-----------------

The argument selector is a special templated metafunction that can return its nth template parameter. This might seem a trivial and unuseful task, however, you should remind yourself that once a template is instantiated there is no way to obtain the actual value of the template parameter from outside the given template. The implementation of the argument selector looks like:

\snippet metaprogramming.hpp ArgumentSelector

In the above code the key role of recursion is clearly followable. The general case of the `select_argument` template class is defined by means of a simple recursion rule, whereas the recursion tail is a partial specialisation of the general case.

Metafunction classes {#metafunclass}
--------------------

Metafunction classes (a.k.a. metafunctors) are classes with a templated inner class that is a metafunction. Just as function classes (a.k.a. functors) have the evaluation operator function that can return a value, metafunction classes have a nested metafunction, that can return a type, as described [above](#metafun). In NiHu, the nested metafunction is always called `apply`.

\snippet metaprogramming.hpp MetaFunctionClass

Placeholders {#placeholder}
------------

Placeholders are special [metafunction classes](#metafunclass) that contain an [argument selector](#argselector) as the nested class. Placeholders do not have any other specialities what would make them differ from other metafunctors, yet, implemenetation of template metaprogramming control sequences and algorithms involve the extensive usage of placeholders, which makes them worth explaining separately.

\snippet metaprogramming.hpp PlaceHolder

For the sake of convenience placeholders are referred to through simplified type definitions.

\snippet metaprogramming.hpp PlaceHolderTypedef

Placeholder expressions {#placeholderexpression}
-----------------------

Placeholder expression is a term referring either to 
	1. a [placeholder](#placeholder)
	2. a template specialisation with at least one argument that is a [placeholder expression](#placeholderexpression)

Deciding whether an expression is a placeholder expression or not, is determined recursively. The snippet below demonstrates how the decision is implemented using variadic templates. (Note: the metafunction called `containsPlaceholderExpression`, which is not demonstrated here, can tell if an argument list contains a placeholder expression.)

\snippet metaprogramming.hpp PlaceHolderDecide

Noticing that the variadic template denoted as `MF<Args...>` refers to a metafunction with the argument `Args`, it is straightforward that the above code implements the decision based on the criteria listed above.

Lambda expressions {#lambdaexpression}
------------------
	
Is a term referring either to 
	1. a metafunction class 
	2. or (as a special case) a placeholder expression.

Lambda metafunctions {#lambdametafun}
--------------------

Are special metafunctions, which take a lambda expression as a template parameter. The purpose of using lambda metafunctions is to turn placeholder expressions into metafunction classes. Lambda metafunctions contain a nested metafunction class.

The apply metafunction {#applymetafun}
----------------------

Is a special metafunction, which invokes the result of a lambda metafunction. In other words, the `apply` metafunction is a shorthand notation.

\snippet metaprogramming.hpp ApplyMetaFunction

Common techniques {#commontechniques}
=================

Type definition forwarding {#typedeffwd}
--------------------------

Type definition forwarding can easily be achieved using inheritance. This is usually applied as a shortcut in order to avoid repetition of typedef lines. One can also make use of the fact that different template specialisations can be derived from different base classes.
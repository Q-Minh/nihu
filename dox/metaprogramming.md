Metaprogramming concepts in NiHu {#metaprogramming}
================================

[TOC]

Introduction {#intro}
============

This document briefly summarises the metaprogramming concepts used in NiHu. 

Nomenclature {#nomen}
============

In the following some metaprogramming elements are explained.

Metafunctions {#metafun}
-------------

Metafunctions are program elements (classes in C++) that can return computed types. Whereas a function returns a value a metafunction returns a type. Metafunctions are implemented in NiHu usings simple structures with nested type definitions, called `type`.

\snippet metaprogramming.hpp MetaFunction

Argument selector {#argselector}
-----------------

The argument selector is a special templated metafunction that can return its nth template parameter.

\snippet metaprogramming.hpp ArgumentSelector

Metafunction classes {#metafunclass}
--------------------

Metafunction classes (a.k.a. metafunctors) are classes with a templated inner class that is a metafunction. In NiHu, the nested metafunction is called `apply`.

\snippet metaprogramming.hpp MetaFunctionClass

Placeholders {#placeholder}
------------

Placeholders are special metafunction classes that contain an argument selector as the nested class.

\snippet metaprogramming.hpp PlaceHolder

For the sake of convenience placeholders are referred to through simplified type definitions.

\snippet metaprogramming.hpp PlaceHolderTypedef

Placeholder expressions {#placeholderexpression}
-----------------------

Placeholder expression is a term referring either to 
	1. a placeholder 
	2. a template specialization with at least one argument that is a placeholder expression

Lambda expressions {#lambdaexpression}
------------------
	
Is a term referring either to 
	1. a metafunction class 
	2. or (as a special case) a placeholder expression.

Lambda metafunctions {#lambdametafun}
--------------------

Is a special metafunction, which takes a lambda expression as a template parameter. The purpose of using lambda metafunctions is to turn placeholder expressions into metafunction classes. Lambda metafunctions contain a nested metafunction class.

The `apply` metafunction {#applymetafun}
------------------------

Is a special metafunction, which invokes the result of a lambda

\snippet metaprogramming.hpp ApplyMetaFunction

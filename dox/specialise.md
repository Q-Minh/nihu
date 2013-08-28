How to specialise an interface class in NiHu {#specialise}
============================================

\page specialise

[CRTP]:http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
[Static polymorphism]:http://en.wikipedia.org/wiki/Template_metaprogramming#Static_polymorphism
[Traits classes]:http://www.cantrip.org/traits.html

[TOC]

Introduction {#specialise_intro}
============

This tutorial explains how to specialise a general interface class in the NiHu toolbox.

All customisable entities in NiHu are implemented using [Static polymorphism], [CRTP] and [Traits Classes].
The main concept of CRTP with traits classes can be summarised as follows:
- The customised class is derived from an interface base class.
- The base class takes the customised class as template argument.
- The base class obtains properties of the custom class from the custom class's traits class, e.g. a traits class template specialised to the custom class.

A minimal example {#specialise_example}
=================

Below you can see a minimal working example where class `entity` is customised.
The naming conventions in NiHu are the following.
If the class to implement would be called `entity`, then the base interface is termed `entity_base`, and the traits class is called `entity_traits`.

The base implementation {#specialise_base}
-----------------------

The code snippet below shows the implementation of the base interface.
This is the part that is already there in the source code.

\snippet specialise.hpp Internal

Carefully reading the code above we can catch all main properties of static polymorphism:
- The base class knows the customised class (although it does not exist yet) from its template parameter `CustomEntity`.
- The base class looks up the traits of the customised class to define the scalar type `scalar_t`.
- The base class defines two interface functions, both of them use the type `scalar_t` in their declarations.
Apparently, the implementation is left for the derived custom class.
- The base class implements an abstract algorithm in terms of the interface functions.
- Functionality implemented in the base class is available in all derived classes, similar to dynamic polymorphism.

Now we are going to customise class `entity_base`. The customised class will be called `custom_entity`.
Customisation is done in three steps.
- First we declare our new class with a forward declaration.
- Then we define the traits class of `custom_entity`
- Finally we define the class `custom_entity` itself

Declaring the custom class {#specialise_declare}
--------------------------

\snippet specialise.hpp Declare

Defining the traits class {#specialise_traits}
-------------------------

The traits class can be defined by explicit template specialization shown below.

\snippet specialise.hpp Traits

We mention that there is no general convenient way to find which fields of the traits class need to be defined.
One possible method is to browse the base class and search for expressions of the form `entity_traits<CustomEntity>::***`.
Other method is to browse the documentation of `entity_base`.
The third method is to read compiler errors.

Defining the custom class {#specialise_class}
-------------------------

The custom class needs to be derived from `entity_base<custom_entity>`.
This relation ensures that the custom class can use all necessary typedefs and functions of the base.
Implementation of the required functions can be done as follows:

\snippet specialise.hpp Class

We are done. After compilation, the abstract algorithm runs with our custom types and implementations.

As before, we mention here that there is no general convenient way to find which interface functions of the base class need to be implemented.
The best practice is to search the documentation and the source code.


Customising with a template {#specialise_template}
===========================

In some cases it is convenient to implement a derived custom class that is a template itself.
There is one minor difference from the above example in this case.
If the derived class is a template, then the type names inherited from the base interface should be fully qualified.
This means that all occurrences of `scalar_t` in the derived class should be replaced by `typename entity_base<custom_entity<T> >::scalar_t`, and the same holds for `scalar_ptr`. One convenient workaround is shown below, retypedefing the definitions of base:

\snippet specialise.hpp Custom template



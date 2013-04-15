Homogeneous collections {#homocollect}
=======================

[TOC]

Polymorphism in NiHu {#nihupoly}
====================

Evaluating the Weighted Residual approach means integrating a kernel over a set of elements extended with some shape functions. Obviously, when programming a BEM, we want to write our integration routine generally, so that our code remains capable to handle as many different element and kernel types as possible.

Dynamic polymorphism {#dynpoly}
--------------------

A straightforward C++ solution to this problem is *dynamic polymorphism* implemented with virtual functions, abstract base classes and heterogeneous collections. Using this technique, an integration routine can be written that receives pointers to abstract kernel and element interfaces classes as input, and performs integration by invoking the functionalities of the abstract interfaces. The specific behaviour of different subproblems are implemented in specialised derived classes of the abstract base.

The main advantage of dynamic polymorphism is that if we build a vector of element pointers, integration over the whole mesh can be performed with a single loop over our container.

The main drawbacks of the dynamic polymorphism approach are the following:
- The actual implementaion functions are dispatched during runtime, so they cannot be inlined by the compiler.
- The size of the result and temporary vectors is unknown compile time, so they need to be stored dynamically, and even simple traversing loops cannot be unrolled by the compiler.

These drawbacks can add up to a significant decrease of overall performance.


Static polymorphism {#statpoly}
--------------------

Static polymorphism means that we implement our integration routine generally

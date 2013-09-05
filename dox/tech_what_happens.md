What happens when NiHu evaluates a weithed residual? {#tech_what_happens}
====================================================

\page tech_what_happens

[TOC]

Introduction {#tech_what_happens_intro}
============

This page explains what happens when one evaluates a weighted double integral

\f$
\displaystyle W_{ij} = \left< w_i, \left(\mathcal{K} d_j\right)_S \right>_F
\f$

or, as coded in a NiHu equaivalent form
~~~~~~~~~~~
W << w * K[d];
~~~~~~~~~~~
where
- `w` and `d` are function spaces created either using function ::create_function_space or using a view generating function like ::constant_view, ::isoparametric_view or ::dirac.
- and `K` is an integral operator created using function ::create_integral_operator from a kernel function.

Operators and proxy objects {#tech_what_happens_operators}
---------------------------

When evaluating the examplified expression, the integral transform (`[]` "index") operator ::integral_operator::operator[] is evaluated first. This function simply returns a temporary ::projection instance.

A projection is a proxy class that simply stores references to the integral operator `K` and the function space `d`. It does not provide any additional data to the integral transform, only functionalities. For example an integral transform can be tested by a test function space, as shown below.

The next step is the evaluation of the testing operator (`*` "multiplication") with a left hand side function space `w` and a right hand side projection (temporary). This testing operation returns a new instance of class ::weighted_residual. This class is again a proxy class storing references to the left hand side function space and the temporary integral transform. The additional functionality provided by the weighted residual class is that it can be evaluated into a matrix.

The next step is calling ::operator<< with the matrix `W` as left hand side and the temporary weighted residual proxy as right hand side. ::operator<< calls the member function template ::weighted_residual::eval of the temporary weighted residual instance. This function is templated on the matrix type, providing that weighted residuals can be evaluated into any indexable type (Eigen matrix of Matlab matrix, for example). This eval function simply forwards the result matrix and the stored test function space reference to the temporary integral transform's member function ::projection::test_on_into parametrised with the test function space `w` and the output matrix reference. This member function _tests_ the integral transform _on_ the test function space and writes the result _into_ the matrix.

Function ::projection::test_on_into is templated on the test function space's type. Therefore, this function is able to determine during compilation if the stored trial function space and the test space parameter are defined on the same mesh type. If yes, then singular integration may be needed. If no, then only regular integrals need to be evaluated. If the mesh types are the same, then the function checks whether the two mesh instances are the same by comparing the mesh objects' addresses.
\note This is why it is important to store function space views in references. If we make copies of meshes, then necessary singular integrals may be omitted during evaluation.

After this decision has been made, the appropriate instance of class template ::assembly gets the scope. This class is responsible for assembling the output matrix from element submatrices, and it is templated on the test function space type, the trial function space type and on the compile time information whether singular integrals need to be computed or not. The static member function ::assembly::eval_into is called with the integral operator, the test function space, the trial function space and the result matrix as arguments.

Let's take a break and summarise what happened up till now. We evaluated some operators, created proxy objects that were used to finally call ::assembly::eval_into parametrised with the two function spaces `w` and `d`, the operator `K` and the result matrix `W`. The only runtime information processed was the comparison of two pointers (if the two mesh types are the same) in order to determine if we need to perform singular integrals. All remaining function calls and proxy instantiations are "optimised out" by the optimising compiler.

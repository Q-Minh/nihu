What happens when NiHu evaluates a weighted residual? {#tech_what_happens}
=====================================================

\page tech_what_happens

[TOC]

Introduction {#tech_what_happens_intro}
============

This page explains what happens when NiHu evaluates a weighted double integral

\f$
\displaystyle W_{ij} = \left< w_i, \left(\mathcal{K} d_j\right)_S \right>_F
\f$

or, as coded in a NiHu equivalent form
~~~~~~~~~~~
W << w * K[d];
~~~~~~~~~~~
where
- `w` and `d` are function spaces created either using function ::create_function_space or using a view generating function like ::constant_view, ::isoparametric_view or ::dirac.
- `K` is an integral operator created using function ::create_integral_operator from a kernel function.
- and `W` is a matrix object

Operators and proxy objects {#tech_what_happens_operators}
===========================

When evaluating the exemplified expression, the integral transform (`[]` "index") operator ::integral_operator::operator[] is evaluated first with the test function space as argument. This function simply returns a temporary ::integral_transform instance. An ::integral_transform is a proxy class that simply stores references to the ::integral_operator `K` and the function_space `d`. The proxy class provides new functionalities: For example an ::integral_transform can be tested by a test ::function_space, as shown below.

The next step is the evaluation of the testing operator (`*` "multiplication") with a left hand side ::function_space `w` and a right hand side ::integral_transform (temporary). This testing operation returns a new proxy object of type ::weighted_residual that stores references to the left hand side ::function_space and the temporary ::integral_transform. The additional functionality provided by the ::weighted_residual class is that it can be evaluated into a matrix.

The next step is calling ::operator<< with the matrix `W` as left hand side and the temporary ::weighted_residual proxy as right hand side. ::operator<< calls the member function template ::weighted_residual::eval with the matrix as argument. This function is templated on the matrix type, so that weighted residuals can be evaluated into any indexable type (Eigen matrix of Matlab matrix, for example). The ::weighted_residual::eval function simply forwards the result matrix reference and the stored test ::function_space reference to the temporary ::integral_transform's member function ::integral_transform::test_on_into that _tests_ the integral transform `K[d]` _on_ the test function space `w` and writes the result _into_ the matrix `W`.

Function ::integral_transform::test_on_into is able to determine during compilation if the stored trial function space `d` and the test function space parameter `w` are defined on the same mesh _type_. If no, then the underlying meshes are surely different and only regular integrals need to be evaluated. If yes, then singular integration may be needed. In the latter case, the mesh instances extracted from the function spaces are compared by address in order to make the decision.
\note This is why it is important to store function space views in references. If we make copies of meshes, then necessary singular integrals may be omitted during evaluation.

After this decision has been made, the appropriate instance of class template ::assembly gets the scope. This class is responsible for assembling the output matrix `W` from element submatrices. Class ::assembly is templated on the test ::function_space type, the trial ::function_space type and on the compile time boolean information whether singular integrals need to be computed or not. The static member function ::assembly::eval_into is called with the ::integral_operator `K`, the test function space `w`, the trial function space `d` and the result matrix `W` as arguments.

Let's take a break and summarise what happened up till now. We evaluated some operators, created proxy objects that were used to finally call ::assembly::eval_into parametrised with the two function spaces `w` and `d`, the operator `K` and the result matrix `W`. The only runtime information processed was the comparison of two pointers (if the two mesh types are the same) in order to determine if we need to perform singular integrals or not. All remaining function calls and proxy instantiations are "optimised out" by the optimising compiler.

Evaluation on homogeneous subspaces {#tech_what_happens_dcall_each}
===================================

Function ::assembly::eval_into splits the test and trial function spaces to homogeneous subspaces and matrix assembly is performed consecutively for each homogeneous subspace pair. For example: if the test function space consists of constant tria and constant quad fields and the trial function space consists of constant and isoparametric quads, then matrix assembly is performed separately for each row of the table below

Test space    | Trial space
--------------|--------------
constant tria | constant quad
constant tria | isoparametric quad
constant quad | constant quad
constant quad | isoparametric quad

This decomposition is performed during compile time, by utilising the library's code generation routine tmp::d_call_each. (d refers to a double nested loop over types)

The assembler's behaviour is slightly different if the integral operator `K` is a local operator (meaning that it acts only within an element, like the identity operator). In this case only those function space pairs are considered where the test and trial function spaces contain the same element type:

Test space    | Trial space
--------------|--------------
constant quad | constant quad
constant quad | isoparametric quad

Now the homogeneous assembler evaluates blocks of the result matrix `W` for each row of the tables, by traversing both homogeneous function spaces and evaluating the double integral

\f$
\displaystyle
W_{ij} = \left< w_i, \left(\mathcal{K}d_j\right)_{S_i}\right>_{F_j}
\f$

where \f$ S_i \f$ and \f$ F_j \f$ denote the subelements of the two (homogeneous sub)domains.


_to be continued_


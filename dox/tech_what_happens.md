What happens when NiHu evaluates a weithed residual? {#tech_what_happens}
====================================================

\page tech_what_happens

[TOC]

Introduction {#tech_what_happens_intro}
============

This page explains what happens when NiHu evaluates a weighted double integral

\f$
\displaystyle W_{ij} = \left< w_i, \left(\mathcal{K} d_j\right)_S \right>_F
\f$

or, as coded in a NiHu equaivalent form
~~~~~~~~~~~
W << w * K[d];
~~~~~~~~~~~
where
- `w` and `d` are function spaces created either using function ::create_function_space or using a view generating function like ::constant_view, ::isoparametric_view or ::dirac.
- `K` is an integral operator created using function ::create_integral_operator from a kernel function.
- and `W` is a matrix object

Operators and proxy objects {#tech_what_happens_operators}
---------------------------

When evaluating the examplified expression, the integral transform (`[]` "index") operator ::integral_operator::operator[] is evaluated first with the test function space as argument. This function simply returns a temporary ::integral_transform instance. An ::integral_transform is a proxy class that simply stores references to the ::integral_operator `K` and the function_space `d`. The proxy class provides new functionalities: For example an ::integral_transform can be tested by a test ::function_space, as shown below.

The next step is the evaluation of the testing operator (`*` "multiplication") with a left hand side ::function_space `w` and a right hand side ::integral_transform (temporary). This testing operation returns a new proxy object of type ::weighted_residual that stores references to the left hand side ::function_space and the temporary ::integral_transform. The additional functionality provided by the ::weighted_residual class is that it can be evaluated into a matrix.

The next step is calling ::operator<< with the matrix `W` as left hand side and the temporary ::weighted_residual proxy as right hand side. ::operator<< calls the member function template ::weighted_residual::eval. This function is templated on the matrix type, so that weighted residuals can be evaluated into any indexable type (Eigen matrix of Matlab matrix, for example). The ::weighted_residual::eval function simply forwards the result matrix reference and the stored test ::function_space reference to the temporary ::integral_transform's member function ::integral_transform::test_on_into that _tests_ the integral transform _on_ the test function space and writes the result _into_ the matrix.

Function ::integral_transform::test_on_into is able to determine during compilation if the stored trial function space and the test function space parameters are defined on the same mesh _type_. If no, then only regular integrals need to be evaluated. If yes, then singular integration may be needed. In the latter case, the mesh instances are compared by address in order to make the decision.
\note This is why it is important to store function space views in references. If we make copies of meshes, then necessary singular integrals may be omitted during evaluation.

After this decision has been made, the appropriate instance of class template ::assembly gets the scope. This class is responsible for assembling the output matrix from element submatrices, and it is templated on the test function space type, the trial function space type and on the compile time information whether singular integrals need to be computed or not. The static member function ::assembly::eval_into is called with the ::integral_operator, the test function space, the trial function space and the result matrix as arguments.

Let's take a break and summarise what happened up till now. We evaluated some operators, created proxy objects that were used to finally call ::assembly::eval_into parametrised with the two function spaces `w` and `d`, the operator `K` and the result matrix `W`. The only runtime information processed was the comparison of two pointers (if the two mesh types are the same) in order to determine if we need to perform singular integrals or not. All remaining function calls and proxy instantiations are "optimised out" by the optimising compiler.

Evaluation on homogeneous subspaces {#tech_what_happens_dcall_each}
-----------------------------------

Function ::assembly::eval_into splits the test and trial function spaces to homogeneous subspaces and matrix assembly is performed for each homogeneous subspace pair. For example: if the test and trial function spaces both consist of constant tria and quad elements, then matrix assembly is performed separately for
- constant tria-constant tria,
- constant tria-constant quad,
- constant quad-constant tria and
- constant quad-constant quad

function space pairs.
This decomposition is performed during compile time, by utilising the library's code generation routine tmp::d_call_each.

The homogeneous assembler makes a new compile time decision. If the integral operator is a local operator (meaning that it acts only within an element, like the identity operator) then processing is only performed if the homogeneous test and trial function spaces contain the same element types:
- constant tria-constant tria
- constant quad-constant quad

In this case, the double code generating loop tmp::d_call_each simplifies automatically to a single (diagonal) loop in compile time.

Now the homogeneous assembler evaluates the result matrix by traversing both function spaces and evaluating the double integral

\f$
\displaystyle
W_{ij} = \left< w_i, \left(\mathcal{K}d_j\right)_{S_i}\right>_{F_j}
\f$

where \f$ S_i \f$ and \f$ F_j \f$ denote the subelements of the two (homogeneous sub)domains.
